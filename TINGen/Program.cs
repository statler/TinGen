using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;

using BitMiracle.LibTiff.Classic;
using MIConvexHull;
using ProjNet.CoordinateSystems;
using ProjNet.CoordinateSystems.Transformations;

internal static class Program
{
    // --------------------------------------------------------------------
    // WCS endpoint (use the OnlineResource href from your GetCapabilities)
    // Your snippet shows staging host; swap to production if/when needed.
    // --------------------------------------------------------------------
    private const string WcsServerUrl =
        "https://services.ga.gov.au/gis/services/DEM_SRTM_1Second_2024/MapServer/WCSServer";

    // Coverage name from your capabilities snippet: <CoverageOfferingBrief><name>1</name>
    private const string CoverageName = "1";

    // Outputs
    private const string OutputGeoTiff = "dem_wcs.tif";
    private const string OutputLandXml = "tin_1km_landxml.xml";

    // Center point (WGS84). Replace with Archer River bridge coordinate.
    private const double CenterLat = -13.4379;
    private const double CenterLon = 142.9416;

    // Area of interest
    private const double RadiusMeters = 1000.0;

    // Controls density: WCS will resample to WIDTH x HEIGHT
    // 1024 gives a dense mesh (can be heavy). 512 is often enough.
    private const int WcsWidth = 1024;
    private const int WcsHeight = 1024;

    private const int MaxTinPoints = 250_000;

    // --------------------------------------------------------------------
    // Projection WKTs (keep multiple and switch TargetProjectionWkt)
    // --------------------------------------------------------------------
    private const string Wkt_GDA94_MGA54 = @"
PROJCS[""GDA94 / MGA zone 54"",
    GEOGCS[""GDA94"",
        DATUM[""Geocentric_Datum_of_Australia_1994"",
            SPHEROID[""GRS 1980"",6378137,298.257222101]],
        PRIMEM[""Greenwich"",0],
        UNIT[""degree"",0.0174532925199433]],
    PROJECTION[""Transverse_Mercator""],
    PARAMETER[""latitude_of_origin"",0],
    PARAMETER[""central_meridian"",141],
    PARAMETER[""scale_factor"",0.9996],
    PARAMETER[""false_easting"",500000],
    PARAMETER[""false_northing"",10000000],
    UNIT[""metre"",1]
]";

    private const string Wkt_GDA2020_MGA56 = @"
PROJCS[""GDA2020 / MGA zone 56"",
    GEOGCS[""GDA2020"",
        DATUM[""Geocentric_Datum_of_Australia_2020"",
            SPHEROID[""GRS 1980"",6378137,298.257222101]],
        PRIMEM[""Greenwich"",0],
        UNIT[""degree"",0.0174532925199433]],
    PROJECTION[""Transverse_Mercator""],
    PARAMETER[""latitude_of_origin"",0],
    PARAMETER[""central_meridian"",153],
    PARAMETER[""scale_factor"",0.9996],
    PARAMETER[""false_easting"",500000],
    PARAMETER[""false_northing"",10000000],
    UNIT[""metre"",1]
]";

    // Select target projection here:
    private static readonly string TargetProjectionWkt = Wkt_GDA94_MGA54;
    // private static readonly string TargetProjectionWkt = Wkt_GDA2020_MGA56;

    public static async Task Main()
    {
        CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;

        // Build coordinate transforms WGS84 <-> target projected CRS
        var ctFactory = new CoordinateTransformationFactory();
        var csFactory = new CoordinateSystemFactory();

        var wgs84 = GeographicCoordinateSystem.WGS84;
        var target = csFactory.CreateFromWkt(TargetProjectionWkt);

        var wgs84ToTarget = ctFactory.CreateFromCoordinateSystems(wgs84, target);
        var targetToWgs84 = ctFactory.CreateFromCoordinateSystems(target, wgs84);

        // Center in target CRS (metres)
        var centerTarget = wgs84ToTarget.MathTransform.Transform(new[] { CenterLon, CenterLat });
        double cx = centerTarget[0];
        double cy = centerTarget[1];

        // 1km square bbox in target CRS, then convert to WGS84 bbox for WCS BBOX
        double minX = cx - RadiusMeters;
        double maxX = cx + RadiusMeters;
        double minY = cy - RadiusMeters;
        double maxY = cy + RadiusMeters;

        var ll = targetToWgs84.MathTransform.Transform(new[] { minX, minY });
        var ur = targetToWgs84.MathTransform.Transform(new[] { maxX, maxY });

        double minLon = ll[0], minLat = ll[1];
        double maxLon = ur[0], maxLat = ur[1];

        Debug.WriteLine($"Center (WGS84): lat={CenterLat}, lon={CenterLon}");
        Debug.WriteLine($"Center (TARGET): X={cx:F3}, Y={cy:F3}");
        Debug.WriteLine($"BBox (WGS84): lon[{minLon},{maxLon}] lat[{minLat},{maxLat}]");
        Debug.WriteLine($"WCS size: {WcsWidth}x{WcsHeight}");

        using var http = CreateHttpClient();

        // Download WCS coverage as GeoTIFF (real elevation values)
        var geotiffBytes = await DownloadWcsGeoTiffAsync(http, minLon, minLat, maxLon, maxLat, WcsWidth, WcsHeight);
        await File.WriteAllBytesAsync(OutputGeoTiff, geotiffBytes);
        Debug.WriteLine($"Downloaded WCS GeoTIFF: {geotiffBytes.Length / 1024.0 / 1024.0:F2} MB -> {OutputGeoTiff}");

        // Read pixels GDAL-free and convert to TIN points in TARGET CRS
        var points = ReadPointsFromGeoTiffBytes(
            geotiffBytes,
            // IMPORTANT: WCS coverage in your capabilities uses CRS84 (lon/lat order).
            // We request CRS=EPSG:4326 but we will map pixels using the bbox we sent.
            minLon, minLat, maxLon, maxLat,
            wgs84ToTarget,
            cx, cy,
            RadiusMeters,
            MaxTinPoints);

        Debug.WriteLine($"TIN points: {points.Count:N0}");

        // Triangulate
        var tri = Triangulate(points);
        Debug.WriteLine($"TIN triangles: {tri.Cells.Count():N0}");

        // Write LandXML in TARGET CRS
        WriteLandXml(OutputLandXml, points, tri);
        Debug.WriteLine($"Wrote LandXML: {Path.GetFullPath(OutputLandXml)}");
    }

    // --------------------------------------------------------------------
    // HTTP helper (prevents GA 403s when UA is missing)
    // --------------------------------------------------------------------
    private static HttpClient CreateHttpClient()
    {
        var http = new HttpClient { Timeout = TimeSpan.FromMinutes(10) };
        http.DefaultRequestHeaders.TryAddWithoutValidation(
            "User-Agent",
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36");
        http.DefaultRequestHeaders.TryAddWithoutValidation("Accept", "image/tiff,application/octet-stream,*/*");
        return http;
    }

    // --------------------------------------------------------------------
    // WCS download (GeoTIFF)
    // Uses WCS 1.0.0 as shown in your capabilities.
    // --------------------------------------------------------------------
    private static async Task<byte[]> DownloadWcsGeoTiffAsync(
        HttpClient http,
        double minLon, double minLat, double maxLon, double maxLat,
        int width, int height)
    {
        // Some ArcGIS WCS instances accept FORMAT=GeoTIFF, others image/tiff.
        // We'll try GeoTIFF first, and if it fails, fall back to image/tiff.
        var bytes = await TryGetCoverage(http, minLon, minLat, maxLon, maxLat, width, height, "GeoTIFF");
        if (LooksLikeTiff(bytes)) return bytes;

        bytes = await TryGetCoverage(http, minLon, minLat, maxLon, maxLat, width, height, "image/tiff");
        if (LooksLikeTiff(bytes)) return bytes;

        string head = string.Join(" ", bytes.Take(16).Select(b => b.ToString("X2")));
        throw new InvalidOperationException($"WCS GetCoverage did not return a TIFF. First bytes={head}");
    }

    private static async Task<byte[]> TryGetCoverage(
        HttpClient http,
        double minLon, double minLat, double maxLon, double maxLat,
        int width, int height,
        string format)
    {
        // NOTE: WCS 1.0.0 uses lon,lat order for CRS84 and often for EPSG:4326 in ArcGIS.
        // Your capabilities lonLatEnvelope uses CRS84. We'll keep the bbox order minLon,minLat,maxLon,maxLat.
        string bbox = $"{minLon.ToString("R")},{minLat.ToString("R")},{maxLon.ToString("R")},{maxLat.ToString("R")}";

        var url =
            $"{WcsServerUrl}?" +
            $"SERVICE=WCS&REQUEST=GetCoverage&VERSION=1.0.0" +
            $"&COVERAGE={Uri.EscapeDataString(CoverageName)}" +
            $"&CRS=EPSG:4326" +
            $"&BBOX={Uri.EscapeDataString(bbox)}" +
            $"&WIDTH={width}&HEIGHT={height}" +
            $"&FORMAT={Uri.EscapeDataString(format)}";

        using var resp = await http.GetAsync(url);
        var bytes = await resp.Content.ReadAsByteArrayAsync();

        if (!resp.IsSuccessStatusCode)
        {
            var snippet = Encoding.UTF8.GetString(bytes.Take(1200).ToArray());
            throw new HttpRequestException($"WCS GetCoverage failed: {(int)resp.StatusCode} {resp.ReasonPhrase}\n{snippet}");
        }

        return bytes;
    }

    // --------------------------------------------------------------------
    // GDAL-free raster read (BitMiracle.LibTiff.NET)
    // Reads single-band Float32 or Int16/UInt16 and converts to float Z.
    //
    // Georeferencing:
    // - We DON'T rely on GeoTIFF tags.
    // - We use the BBOX we requested + WIDTH/HEIGHT to map pixel centers to lon/lat.
    // --------------------------------------------------------------------
    private static List<TinPoint> ReadPointsFromGeoTiffBytes(
    byte[] tiffBytes,
    double minLon, double minLat, double maxLon, double maxLat,
    ICoordinateTransformation wgs84ToTarget,
    double cx, double cy,
    double radiusMeters,
    int maxPoints)
    {
        double r2 = radiusMeters * radiusMeters;

        using var ms = new MemoryStream(tiffBytes);
        ms.Position = 0;
        using var tif = Tiff.ClientOpen("mem", "r", ms, new TiffStream());

        if (tif == null) throw new InvalidOperationException("Could not open TIFF (ClientOpen returned null).");

        int width = GetIntTag(tif, TiffTag.IMAGEWIDTH, 0);
        int height = GetIntTag(tif, TiffTag.IMAGELENGTH, 0);
        int bitsPerSample = GetIntTag(tif, TiffTag.BITSPERSAMPLE, 0);
        int sampleFormat = GetIntTag(tif, TiffTag.SAMPLEFORMAT, (int)SampleFormat.IEEEFP);
        int samplesPerPixel = GetIntTag(tif, TiffTag.SAMPLESPERPIXEL, 1);

        Debug.WriteLine($"TIFF tiled={tif.IsTiled()} width={width} height={height}");
        Debug.WriteLine($"Bits={bitsPerSample} SampleFormat={sampleFormat} SamplesPerPixel={samplesPerPixel}");

        if (width <= 0 || height <= 0)
            throw new InvalidOperationException($"Invalid TIFF size. width={width}, height={height}");

        if (samplesPerPixel != 1)
            throw new InvalidOperationException($"Expected single-band coverage. SamplesPerPixel={samplesPerPixel}");

        // Geo-referencing from request bbox
        double dLon = (maxLon - minLon) / width;
        double dLat = (maxLat - minLat) / height;

        // Bytes per sample (single band)
        int bytesPerSample = bitsPerSample switch
        {
            8 => 1,
            16 => 2,
            32 => 4,
            _ => throw new InvalidOperationException($"Unsupported BitsPerSample={bitsPerSample}")
        };

        var pts = new List<TinPoint>(Math.Min(width * height, maxPoints));

        if (tif.IsTiled())
        {
            int tileWidth = GetIntTag(tif, TiffTag.TILEWIDTH, 0);
            int tileHeight = GetIntTag(tif, TiffTag.TILELENGTH, 0);
            if (tileWidth <= 0 || tileHeight <= 0)
                throw new InvalidOperationException($"TIFF says tiled but missing tile dims. tileWidth={tileWidth}, tileHeight={tileHeight}");

            int tileSize = tif.TileSize();
            var tileBuf = new byte[tileSize];

            for (int ty = 0; ty < height; ty += tileHeight)
            {
                for (int tx = 0; tx < width; tx += tileWidth)
                {
                    Array.Clear(tileBuf, 0, tileBuf.Length);

                    // Read tile (band/sample fixed at 0 for single-band)
                    int read = tif.ReadEncodedTile(tif.ComputeTile(tx, ty, 0, 0), tileBuf, 0, tileSize);
                    if (read <= 0)
                        throw new InvalidOperationException($"Failed reading tile at ({tx},{ty}).");

                    int rowsInTile = Math.Min(tileHeight, height - ty);
                    int colsInTile = Math.Min(tileWidth, width - tx);

                    // Tile buffer layout is row-major for the tile
                    for (int rowInTile = 0; rowInTile < rowsInTile; rowInTile++)
                    {
                        int y = ty + rowInTile;
                        double lat = maxLat - (y + 0.5) * dLat;

                        int rowOffset = rowInTile * tileWidth * bytesPerSample;

                        for (int colInTile = 0; colInTile < colsInTile; colInTile++)
                        {
                            int x = tx + colInTile;
                            double lon = minLon + (x + 0.5) * dLon;

                            int offset = rowOffset + colInTile * bytesPerSample;
                            float z = ReadSampleAsFloat(tileBuf, offset, bitsPerSample, sampleFormat);
                            if (float.IsNaN(z)) continue;

                            var xy = wgs84ToTarget.MathTransform.Transform(new[] { lon, lat });
                            double px = xy[0];
                            double py = xy[1];

                            double dx = px - cx;
                            double dy = py - cy;
                            if (dx * dx + dy * dy <= r2)
                            {
                                pts.Add(new TinPoint(px, py, z));
                                if (pts.Count >= maxPoints) return DedupXY(pts);
                            }
                        }
                    }
                }
            }
        }
        else
        {
            // Stripped TIFF: read strip-by-strip
            int stripSize = tif.StripSize();
            int rowsPerStrip = GetIntTag(tif, TiffTag.ROWSPERSTRIP, 0);
            if (rowsPerStrip <= 0) rowsPerStrip = 1;

            var stripBuf = new byte[stripSize];
            int strips = tif.NumberOfStrips();

            int yBase = 0;

            for (int s = 0; s < strips; s++)
            {
                Array.Clear(stripBuf, 0, stripBuf.Length);

                int read = tif.ReadEncodedStrip(s, stripBuf, 0, stripSize);
                if (read <= 0)
                    throw new InvalidOperationException($"Failed reading strip {s}.");

                int rowsInStrip = Math.Min(rowsPerStrip, height - yBase);

                for (int rowInStrip = 0; rowInStrip < rowsInStrip; rowInStrip++)
                {
                    int y = yBase + rowInStrip;
                    double lat = maxLat - (y + 0.5) * dLat;

                    int rowOffset = rowInStrip * width * bytesPerSample;

                    for (int x = 0; x < width; x++)
                    {
                        double lon = minLon + (x + 0.5) * dLon;

                        int offset = rowOffset + x * bytesPerSample;
                        float z = ReadSampleAsFloat(stripBuf, offset, bitsPerSample, sampleFormat);
                        if (float.IsNaN(z)) continue;

                        var xy = wgs84ToTarget.MathTransform.Transform(new[] { lon, lat });
                        double px = xy[0];
                        double py = xy[1];

                        double dx = px - cx;
                        double dy = py - cy;
                        if (dx * dx + dy * dy <= r2)
                        {
                            pts.Add(new TinPoint(px, py, z));
                            if (pts.Count >= maxPoints) return DedupXY(pts);
                        }
                    }
                }

                yBase += rowsPerStrip;
            }
        }

        return DedupXY(pts);
    }

    private static float ReadSampleAsFloat(byte[] buffer, int offset, int bitsPerSample, int sampleFormat)
    {
        if (bitsPerSample == 32 && sampleFormat == (int)SampleFormat.IEEEFP)
            return BitConverter.ToSingle(buffer, offset);

        if (bitsPerSample == 16 && sampleFormat == (int)SampleFormat.INT)
            return BitConverter.ToInt16(buffer, offset);

        if (bitsPerSample == 16 && sampleFormat == (int)SampleFormat.UINT)
            return BitConverter.ToUInt16(buffer, offset);

        if (bitsPerSample == 8 && (sampleFormat == 0 || sampleFormat == (int)SampleFormat.UINT))
            return buffer[offset];

        throw new InvalidOperationException($"Unsupported pixel type. BitsPerSample={bitsPerSample}, SampleFormat={sampleFormat}");
    }


    private static int GetIntTag(Tiff tif, TiffTag tag, int defaultValue)
    {
        var field = tif.GetField(tag);
        if (field == null || field.Length == 0) return defaultValue;
        return field[0].ToInt();
    }

    private static bool LooksLikeTiff(byte[] bytes)
    {
        if (bytes.Length < 4) return false;
        return (bytes[0] == 0x49 && bytes[1] == 0x49 && bytes[2] == 0x2A && bytes[3] == 0x00) ||
               (bytes[0] == 0x4D && bytes[1] == 0x4D && bytes[2] == 0x00 && bytes[3] == 0x2A);
    }

    private static List<TinPoint> DedupXY(List<TinPoint> points)
    {
        // Quantize to 1 mm for stable hashing in projected metres
        const double q = 0.001;
        var seen = new HashSet<(long, long)>();
        var result = new List<TinPoint>(points.Count);

        foreach (var p in points)
        {
            long kx = (long)Math.Round(p.X / q);
            long ky = (long)Math.Round(p.Y / q);
            if (seen.Add((kx, ky)))
                result.Add(p);
        }

        return result;
    }

    // --------------------------------------------------------------------
    // MIConvexHull triangulation
    // --------------------------------------------------------------------
    private sealed class TinVertex : IVertex
    {
        public double[] Position { get; }  // [X,Y] in target CRS
        public int Index { get; }          // original points index

        public TinVertex(double x, double y, int index)
        {
            Position = new[] { x, y };
            Index = index;
        }
    }

    private sealed class TinCell : TriangulationCell<TinVertex, TinCell> { }

    private static DelaunayTriangulation<TinVertex, TinCell> Triangulate(List<TinPoint> points)
    {
        if (points == null) throw new ArgumentNullException(nameof(points));
        if (points.Count < 3) throw new ArgumentException("Need at least 3 points to triangulate.", nameof(points));

        var verts = new List<TinVertex>(points.Count);
        for (int i = 0; i < points.Count; i++)
            verts.Add(new TinVertex(points[i].X, points[i].Y, i));

        return DelaunayTriangulation<TinVertex, TinCell>.Create(verts,  1e-9);
    }

    // --------------------------------------------------------------------
    // LandXML writer in target projected CRS (X/Y metres, Z elevation)
    // --------------------------------------------------------------------
    private static void WriteLandXml(string path, List<TinPoint> points, DelaunayTriangulation<TinVertex, TinCell> tri)
    {
        using var sw = new StreamWriter(path, false, new UTF8Encoding(false));

        sw.WriteLine(@"<?xml version=""1.0"" encoding=""UTF-8""?>");
        sw.WriteLine(@"<LandXML xmlns=""http://www.landxml.org/schema/LandXML-1.2"" version=""1.2"">");
        sw.WriteLine(@"  <Surfaces>");
        sw.WriteLine($@"    <Surface name=""TIN_1km"" desc=""WCS {CoverageName} DEM_SRTM_1Second_2024 within 1km radius"">");
        sw.WriteLine(@"      <Definition surfType=""TIN"">");
        sw.WriteLine(@"        <Pnts>");

        for (int i = 0; i < points.Count; i++)
        {
            int id = i + 1;
            var p = points[i];
            sw.WriteLine($@"          <P id=""{id}"">{p.X.ToString("F3", CultureInfo.InvariantCulture)} {p.Y.ToString("F3", CultureInfo.InvariantCulture)} {p.Z.ToString("F3", CultureInfo.InvariantCulture)}</P>");
        }

        sw.WriteLine(@"        </Pnts>");
        sw.WriteLine(@"        <Faces>");

        foreach (var cell in tri.Cells)
        {
            int a = cell.Vertices[0].Index + 1;
            int b = cell.Vertices[1].Index + 1;
            int c = cell.Vertices[2].Index + 1;

            if (a == b || b == c || a == c)
                continue;

            sw.WriteLine($@"          <F>{a} {b} {c}</F>");
        }

        sw.WriteLine(@"        </Faces>");
        sw.WriteLine(@"      </Definition>");
        sw.WriteLine(@"    </Surface>");
        sw.WriteLine(@"  </Surfaces>");
        sw.WriteLine(@"</LandXML>");
    }

    private readonly record struct TinPoint(double X, double Y, double Z);
}
