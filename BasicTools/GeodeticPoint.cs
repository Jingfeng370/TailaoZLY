using System;
using System.Collections.Generic;
using System.Text;

namespace TailaoZLY.BasicTools
{
    public class GeodeticPoint
    {
        #region

        private  ReferenceEllipsoid re = ReferenceEllipsoid.CGCS2000;
        public  ReferenceEllipsoid RefEllipsoid { set { re = value; } get { return re; } }

        private double gX = 0;
        private double gY = 0;
        private double gZ = 0;
        public double X { get { return gX; } }
        public double Y { get { return gY; } }
        public double Z { get { return gZ; } }


        private double gB = 0;
        private double gL = 0;
        private double gH = 0;
        public double B { get { return gB; } }
        public double L { get { return gL; } }
        /// <summary>
        ///大地高
        /// </summary>
        public double H { get { return gH; } }

        private double gx = 0;
        private double gy = 0;
        public double x { get { return gx; } }
        public double y { get { return gy; } }

        private double gh = 0;
        /// <summary>
        /// 正常高
        /// </summary>
        public double h { get { return gh; } set { gh = value; } }

        private double gKsi = 0;
        /// <summary>
        /// 高程异常
        /// </summary>
        public double Ksi { get { return gKsi; } set { gKsi = value; } }
        #endregion
        public GeodeticPoint(ReferenceEllipsoid refEll, GeoCoordinateType geoType, params string[] data)
        {
            re = refEll;
            GeoXYZ gXYZ = new GeoXYZ(0, 0, 0);
            GeoBLH gBLH = new GeoBLH(0, 0, 0);
            GeoxyH gxyH = new GeoxyH(0, 0, 0);
            switch (geoType)
            {
                case GeoCoordinateType.XYZ:
                    if (data.Length != 3) break;
                    else
                    {
                        gXYZ = new GeoXYZ(data[0], data[1], data[2]);
                        gBLH = ConvertXYZToBLH(gXYZ);
                        gxyH = ConvertBLHToxyH(gBLH);
                    }
                    break;
                case GeoCoordinateType.BLH:
                    if (data.Length == 2)
                    {
                        gBLH = new GeoBLH(data[0], data[1]);
                        gXYZ = ConvertBLHToXYZ(gBLH);
                        gxyH = ConvertBLHToxyH(gBLH);
                    }
                    else if (data.Length == 3)
                    {
                        gBLH = new GeoBLH(data[0], data[1], data[2]);
                        gXYZ = ConvertBLHToXYZ(gBLH);
                        gxyH = ConvertBLHToxyH(gBLH);
                    }
                    else
                        break;
                    break;
                default:
                    if (data.Length == 2)
                    {
                        gxyH = new GeoxyH(data[0], data[1]);
                        gBLH = ConvertxyHToBLH(gxyH);
                        gXYZ = ConvertBLHToXYZ(gBLH);
                    }
                    else if (data.Length == 3)
                    {
                        gxyH = new GeoxyH(data[0], data[1],data[3]);
                        gBLH = ConvertxyHToBLH(gxyH);
                        gXYZ = ConvertBLHToXYZ(gBLH);
                    }
                    else
                        break;
                    break;
            }
            gX = gXYZ.gX;
            gY = gXYZ.gY;
            gZ = gXYZ.gZ;
            gB = gBLH.gB;
            gL = gBLH.gL;
            gH = gBLH.gH;
            gx = gxyH.gx;
            gy = gxyH.gy;
        }

        /// <summary>
        /// 将大地坐标转换为高斯平面坐标
        /// </summary>
        /// <param name="gBLH">大地坐标</param>
        /// <returns>高斯平面坐标</returns>
        private GeoxyH ConvertBLHToxyH(GeoBLH gBLH)
        {
            double _h = gBLH.gH;
            double X = GeodeticCalculate.LengthOfMeridian(re, gBLH.gB);
            double N = GeodeticCalculate.UnitaryCircleRadius(re, gBLH.gB);
            double centerL = (gBLH.gL + 3600 * 360.0) % (3600 * 360.0);
            int beltNo = (int)(centerL / (6.0 * 3600)) + 1;
            centerL = 6.0 * 3600 * beltNo - 3.0 * 3600;
            double dl = gBLH.gL - centerL;
            double sinB = Math.Sin(gBLH.gB / Angle.Rou);
            double cosB = Math.Cos(gBLH.gB / Angle.Rou);
            double t2 = Math.Pow(Math.Tan(gBLH.gB / Angle.Rou), 2);
            double yita2 = re.firste2 * cosB * cosB;
            double x1 = N * sinB * cosB * dl * dl / (2.0 * Angle.Rou * Angle.Rou);
            double x2 = N * sinB * Math.Pow(cosB, 3) * (5 - t2 + 9 * yita2 + 4 * yita2 * yita2) * Math.Pow(dl, 4) / (24.0 * Math.Pow(Angle.Rou, 4));
            double x3 = N * sinB * Math.Pow(cosB, 5) * (61 - 58 * t2 + t2 * t2) * Math.Pow(dl, 6) / (720.0 * Math.Pow(Angle.Rou, 6));
            double y1 = N * cosB * dl / Angle.Rou;
            double y2 = N * Math.Pow(cosB, 3) * (1 - t2 + yita2) * Math.Pow(dl, 3) / (6.0 * Math.Pow(Angle.Rou, 3));
            double y3 = N * Math.Pow(cosB, 5) * (5 - 18 * t2 + t2 * t2 + 14 * yita2 - 58 * t2 * yita2) * Math.Pow(dl, 5) / (120.0 * Math.Pow(Angle.Rou, 5));
            double _x = X + x1 + x2 + x3;
            double _y = y1 + y2 + y3 + 500000.0 + beltNo * 1000000.0;
            return new GeoxyH(_x, _y, _h);
        }
        /// <summary>
        /// 将空间直角坐标转换为大地坐标
        /// </summary>
        /// <param name="gXYZ">空间直角坐标</param>
        /// <returns>大地坐标</returns>
        private GeoBLH ConvertXYZToBLH(GeoXYZ gXYZ)
        {
            double l = Math.Atan2(gXYZ.gY, gXYZ.gX);
            l = l > 0 ? l : 2 * Math.PI + l;
            double tanB0 = gXYZ.gZ / Math.Sqrt(gXYZ.gX * gXYZ.gX + gXYZ.gY * gXYZ.gY);
            double tanB = 0;
            do
            {
                tanB = (gXYZ.gZ + re.a * re.firste2 * tanB0 / Math.Sqrt(1 + tanB0 * tanB0 - re.firste2 * tanB0 * tanB0)) / Math.Sqrt(gXYZ.gX * gXYZ.gX + gXYZ.gY * gXYZ.gY);
                double temp = tanB0;
                tanB0 = tanB;
                tanB = temp;
            }
            while (Math.Abs(tanB - tanB0) > 1.0e-11);
            double _B = Math.Atan(tanB) * Angle.Rou;
            double _L = l * Angle.Rou;
            double n = GeodeticCalculate.UnitaryCircleRadius(re, _B);
            double _H = Math.Sqrt(gXYZ.gX * gXYZ.gX + gXYZ.gY * gXYZ.gY) / Math.Cos(_B / Angle.Rou) - n;
            return new GeoBLH(_B, _L, _H);
        }
        /// <summary>
        /// 将高斯平面坐标转化为大地坐标
        /// </summary>
        /// <param name="gxyH">高斯平面坐标</param>
        /// <returns>大地坐标</returns>
        private GeoBLH ConvertxyHToBLH(GeoxyH gxyH)
        {
            int beltNo = (int)(gxyH.gy / 1000000.0);
            double y = gxyH.gy % 1000000.0 - 500000.0;
            double Bf = GeodeticCalculate.LengthOfMeridianAntiCaculate(re, gxyH.gx);
            double tf = Math.Tan(Bf / Angle.Rou);
            double cosbf = Math.Cos(Bf / Angle.Rou);
            double yitaf = re.firste2 * cosbf * cosbf;
            double Mf = GeodeticCalculate.MeridianCircleRadius(re, Bf);
            double Nf = GeodeticCalculate.UnitaryCircleRadius(re, Bf);
            double b1 = -(tf * Math.Pow(y, 2) / (2.0 * Mf * Nf));
            double b2 = tf * (5 + 3 * Math.Pow(tf, 2) + yitaf - 9 * yitaf * tf * tf) * Math.Pow(y, 4) / (24 * Mf * Math.Pow(Nf, 3));
            double b3 = -(tf * Math.Pow(y, 6)) * (61 + 90 * tf * tf + 45 * tf * tf * tf * tf) / (720 * Mf * Math.Pow(Nf, 5));
            double B = Bf / Angle.Rou + b1 + b2 + b3;
            double L = y / (Nf * cosbf) - Math.Pow(y, 3) * (1 + 2 * tf * tf + yitaf) / (6 * Math.Pow(Nf, 3) * cosbf) + Math.Pow(y, 5) * (5 + 28 * tf * tf + 24 * Math.Pow(tf, 4) + 6 * yitaf + 8 * yitaf * tf * tf) / (120 * Math.Pow(Nf, 5) * cosbf);
            double _B = Bf + (b1 + b2 + b3) * Angle.Rou;
            double _L = L * Angle.Rou + beltNo * 6.0 * 3600 - 3.0 * 3600;
            double _H = gxyH.gH;
            return new GeoBLH(_B, _L, _H);
        }
        /// <summary>
        /// 将大地坐标转换为空间直角坐标
        /// </summary>
        /// <param name="gblh">大地坐标</param>
        /// <returns>空间直角坐标</returns>
        private GeoXYZ ConvertBLHToXYZ(GeoBLH gblh)
        {
            double n = GeodeticCalculate.UnitaryCircleRadius(re, gblh.gB);
            double cosB = Math.Cos(gblh.gB / Angle.Rou);
            double sinB = Math.Sin(gblh.gB / Angle.Rou);
            double cosL = Math.Cos(gblh.gL / Angle.Rou);
            double sinL = Math.Sin(gblh.gL / Angle.Rou);
            double _X = (n + gblh.gH) * cosB * cosL;
            double _Y = (n + gblh.gH) * cosB * sinL;
            double _Z = (n * (1 - re.firste2) + gblh.gH) * sinB;
            return new GeoXYZ(_X, _Y, _Z);
        }
    

        public GeodeticPoint(ReferenceEllipsoid refEll, GeoCoordinateType geoType, params double[] data)
        {
            re = refEll;
            GeoXYZ gXYZ = new GeoXYZ(0, 0, 0);
            GeoBLH gBLH = new GeoBLH(0, 0, 0);
            GeoxyH gxyH = new GeoxyH(0, 0, 0);
            switch (geoType)
            {
                case GeoCoordinateType.XYZ:
                    if (data.Length != 3) break;
                    else
                    {
                        gXYZ = new GeoXYZ(data[0], data[1], data[2]);
                        gBLH = ConvertXYZToBLH(gXYZ);
                        gxyH = ConvertBLHToxyH(gBLH);
                    }
                    break;
                case GeoCoordinateType.BLH:
                    if (data.Length == 2)
                    {
                        gBLH = new GeoBLH(data[0], data[1]);
                        gXYZ = ConvertBLHToXYZ(gBLH);
                        gxyH = ConvertBLHToxyH(gBLH);
                    }
                    else if (data.Length == 3)
                    {
                        gBLH = new GeoBLH(data[0], data[1], data[2]);
                        gXYZ = ConvertBLHToXYZ(gBLH);
                        gxyH = ConvertBLHToxyH(gBLH);
                    }
                    else
                        break;
                    break;
                default:
                    if (data.Length == 2)
                    {
                        gxyH = new GeoxyH(data[0], data[1]);
                        gBLH = ConvertxyHToBLH(gxyH);
                        gXYZ = ConvertBLHToXYZ(gBLH);
                    }
                    else if (data.Length == 3)
                    {
                        gxyH = new GeoxyH(data[0], data[1], data[3]);
                        gBLH = ConvertxyHToBLH(gxyH);
                        gXYZ = ConvertBLHToXYZ(gBLH);
                    }
                    else
                        break;
                    break;
            }
            gX = gXYZ.gX;
            gY = gXYZ.gY;
            gZ = gXYZ.gZ;
            gB = gBLH.gB;
            gL = gBLH.gL;
            gH = gBLH.gH;
            gx = gxyH.gx;
            gy = gxyH.gy;
        }
        private class GeoXYZ
        {

            public double gX = 0;
            public double gY = 0;
            public double gZ = 0;
            public GeoXYZ(double gx, double gy, double gz)
            {
                gX = gx;
                gY = gy;
                gZ = gz;
            }
            public GeoXYZ(string gx, string gy, string gz)
            {
                double.TryParse(gx, out gX);
                double.TryParse(gy, out gY);
                double.TryParse(gz, out gZ);
            }
        }
        private class GeoBLH
        {
            public double gB = 0;
            public double gL = 0;
            public double gH = 0;
            /// <summary>
            /// 由浮点数生成大地坐标，以秒为单位
            /// </summary>
            /// <param name="b"></param>
            /// <param name="l"></param>
            /// <param name="h"></param>
            public GeoBLH(double b, double l, double h = 0)
            {
                gB = b ;
                gL = l ;
                gH = h;
            }
            public GeoBLH(string b, string l, string h = "0")
            {
                gB = new Angle(b, AngleStyle.Dms).Seconds;
                gL = new Angle(l, AngleStyle.Dms).Seconds;
                double.TryParse(h, out gH);
            }
        }
        private class GeoxyH
        {
            public double gx = 0;
            public double gy = 0;
            public double gH = 0;

            public GeoxyH(double x, double y, double h = 0)
            {
                gx = x;
                gy = y;
                gH = h;
            }
            public GeoxyH(string x, string y, string h = "0")
            {
                double.TryParse(x, out gx);
                double.TryParse(y, out gy);
                double.TryParse(h, out gH);
            }
        }


        
    }
    public enum GeoCoordinateType
    {
        XYZ,
        BLH,
        xyH
    }
}

