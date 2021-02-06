using System;
using System.Collections.Generic;
using System.Text;

namespace TailaoZLY.BasicTools
{
    public static class GeodeticCalculate
    {

        /// <summary>
        /// 文森特大地公式正算
        /// </summary>
        /// <param name="blh">起始点大地坐标</param>
        /// <param name="ang">起始点大地方位角,以秒为单位</param>
        /// <param name="s">大地距离</param>
        /// <param name="alfa2">终点大地方位角A</param>
        /// <returns>终点大地坐标</returns>
        public static GeodeticPoint VincentyCalculate(GeodeticPoint blh, Angle ang, double s, out Angle alfa2)
        {
            //判断 b1,l1,a1取值是否合法
            double b1, l1, a1;
            #region
            b1 = blh.B <= 90.0 * 3600 ? blh.B : 90.0 * 3600;
            b1 = blh.B >= -90.0 * 3600 ? blh.B : -90.0 * 3600;
            l1 = blh.L <= 180.0 * 3600 ? blh.L : 360.0 * 3600 - blh.L;
            l1 = blh.L >= -180.0 * 3600 ? blh.L : 360.0 * 3600 + blh.L;
            a1 = ang.Seconds < 0 ? ang.Seconds + 360.0 * 3600 : ang.Seconds;
            b1 = b1 /Angle.Rou;
            l1 = l1 / Angle.Rou;
            a1 = a1 / Angle.Rou;
            #endregion
            //由起始点计算球面归化纬度
            double u1 = 0;
            if (Math.Abs(b1) == Math.PI * 0.5)
            {
                u1 = b1;
            }
            else
            {
                u1 = Math.Atan((1 - blh.RefEllipsoid.f) * Math.Tan(b1));
            }
            //计算球面角距delta1；
            #region
            double delta1 = 0;
            if (u1 == 0)
            {
                delta1 = 0;
            }
            else if (Math.Abs(u1) == Math.PI / 2.0)
            {
                delta1 = Math.PI * 0.5;
            }
            else if (a1 == Math.PI * 0.5 || a1 == Math.PI * 1.5)
            {
                delta1 = Math.PI * 0.5;
            }
            else
            {
                double x = Math.Atan2(Math.Tan(u1), Math.Cos(a1));
                if (x >= 0)
                {
                    delta1 = x;
                }
                else
                {
                    delta1 = x + Math.PI;
                }
            }
            #endregion
            //计算临界点的球面归化纬度
            #region
            double un = 0;
            double y = Math.Acos(Math.Cos(u1) * Math.Abs(Math.Sin(a1)));
            if (u1 < 0)
            {
                un = -y;
            }
            else if (u1 == 0 && Math.Cos(a1) < 0)
            {
                un = -y;
            }
            else
                un = y;
            #endregion
            // 计算球面角距
            //double t = ep.ee2 * Math.Sin(un) * Math.Sin(un) / 4.0;
            double t = blh.RefEllipsoid.seconde2 * Math.Sin(un) * Math.Sin(un) / 4.0;
            double K1 = 1 + t * (1 - t * (3 - t * (5 - 11 * t)) / 4.0);
            double K2 = t * (1 - t * (2 - t * (37 - 94 * t) / 8.0));
            double delta = 0;
            double deltam = 0;
            double Ddelta1 = 0;
            double Ddelta2 = 0;
            do
            {
                delta = s / (K1 * blh.RefEllipsoid.b) + Ddelta1;
                deltam = 2 * delta1 + delta;
                Ddelta2 = K2 * Math.Sin(delta) * (Math.Cos(deltam) + K2 * (Math.Cos(delta) * Math.Cos(2 * deltam) + K2 * (1 + 2 * Math.Cos(2 * delta)) * Math.Cos(3 * deltam) / 6) / 4);
                double temp = Ddelta2;
                Ddelta2 = Ddelta1;
                Ddelta1 = temp;
            } while (Math.Abs(Ddelta1 - Ddelta2) > 1.0E-13);
            //计算P2点规划纬度u2,然后计算大地纬度B2
            double u2 = Math.Asin(Math.Sin(un) * Math.Sin(delta + delta1));
            double B2 = 0;
            if (Math.Abs(u2) == Math.PI * 0.5)
            {
                B2 = Math.PI * 0.5;
            }
            else
            {
                B2 = Math.Atan2(Math.Tan(u2), 1 - blh.RefEllipsoid.f);
            }
            //计算P2点的大地方位角
            #region
            double A2 = 0;
            if (a1 == 0 || a1 == Math.PI)
            {
                double z = Math.Cos(delta1) * Math.Cos(delta1 + delta);
                if (z < 0)
                {
                    A2 = a1;
                }
                else if (z > 0)
                {
                    A2 = a1 + Math.PI;
                }
                else
                {
                    if (Math.Abs(u1) == 0.5 * Math.PI && u1 == u2)
                    {
                        A2 = a1;
                    }
                    if (Math.Abs(u1) == Math.PI)
                    {
                        A2 = a1 + Math.PI;
                    }
                    if (Math.Abs(u2) == 0.5 * Math.PI && u1 != u2)
                    {
                        A2 = a1 + Math.PI;
                    }
                }
            }
            else
            {
                double m = Math.Sin(a1) * Math.Sin(u1) / Math.Cos(u2);
                double n = (Math.Cos(delta) * Math.Sin(u2) - Math.Sin(u1)) / (Math.Sin(delta) * Math.Cos(u2));
                if (Math.Sin(delta) == 0)
                {
                    if (Math.Cos(delta) == 1)
                    {
                        A2 = a1 + Math.PI;
                    }
                    if (Math.Cos(delta) == -1)
                    {
                        A2 = 2.0 * Math.PI - a1;
                    }
                }
                else
                {
                    if (m >= 0)
                    {
                        A2 = Math.PI + Math.Acos(n);
                    }
                    else
                    {
                        A2 = 3.0 * Math.PI - Math.Acos(n);
                    }
                }

            }
            while (A2 >= 2.0 * Math.PI)
            {
                A2 -= 2.0 * Math.PI;
            }
            while (A2 < 0)
            {
                A2 += 2.0 * Math.PI;
            }
            #endregion
            //计算球面经差和经差改正项
            #region
            double w = 0;
            double dw = 0;
            if (a1 == 0 || a1 == Math.PI)
            {
                double z = Math.Cos(delta1) * Math.Cos(delta1 + delta);
                if (z >= 0)
                    w = 0;
                else
                    w = Math.PI;
            }
            else
            {
                double v = blh.RefEllipsoid.f * Math.Sin(un) * Math.Sin(un) / 4.0;
                double K3 = v * (1 + blh.RefEllipsoid.f + blh.RefEllipsoid.f * blh.RefEllipsoid.f - v * (3 + 7 * blh.RefEllipsoid.f - 13 * v));
                double d = (1 - K3) * blh.RefEllipsoid.f * Math.Cos(un) * (delta + K3 * Math.Sin(delta) * (Math.Cos(deltam) + K3 * Math.Cos(delta) * Math.Cos(2.0 * deltam)));
                double m = Math.Sin(delta) * Math.Sin(a1) / Math.Cos(u2);
                double n = (Math.Cos(delta) - Math.Sin(u1) * Math.Sin(u2)) / (Math.Cos(u1) * Math.Cos(u2));
                if (Math.Sin(delta) == 0)
                {
                    dw = 0;
                    if (Math.Cos(delta) == 1)
                    {
                        w = 0;
                    }
                    else
                        w = Math.PI;
                }
                else
                {
                    if (m >= 0)
                    {
                        dw = d;
                        w = Math.Acos(n);
                    }
                    else
                    {
                        dw = -d;
                        w = 2.0 * Math.PI - Math.Acos(n);
                    }
                }
                while (w > Math.PI)
                    w -= 2.0 * Math.PI;
                while (w < -Math.PI)
                    w += 2.0 * Math.PI;
            }
            #endregion
            //计算大地经度
            double L2 = l1 + w - dw;
            while (L2 > 2.0 * Math.PI)
                L2 -= 2.0 * Math.PI;
            while (L2 < -2.0 * Math.PI)
                L2 += 2.0 * Math.PI;
            alfa2 = new Angle(A2 * Angle.Rou);
            return new GeodeticPoint(blh.RefEllipsoid, GeoCoordinateType.BLH, B2 * Angle.Rou, L2 * Angle.Rou);
        }
        /// <summary>
        /// 由文森特公式进行大地反算
        /// </summary>
        /// <param name="blh1">起点坐标</param>
        /// <param name="blh2">终点坐标</param>
        /// <param name="ang1">正方位角</param>
        /// <param name="ang2">反方位角</param>
        /// <returns>大地线长</returns>
        public static double VincentyAntiCalculate(GeodeticPoint blh1, GeodeticPoint blh2, out Angle ang1, out Angle ang2)
        {
            double[] res = new double[3];
            double A1 = 0;
            double A2 = 0;
            double S = 0;
            double L = 0;
            double B1, L1, B2, L2;
            //对B1，B2，L1，L2规范化处理
            #region
            B1 = blh1.B <= 90.0 * 3600 ? blh1.B : 90.0 * 3600;
            B1 = blh1.B >= -90.0 * 3600 ? blh1.B : -90.0 * 3600;
            L1 = blh1.L <= 180.0 * 3600 ? blh1.L : blh1.L - 360.0 * 3600;
            L1 = blh1.L >= -180.0 * 3600 ? blh1.L : 360.0 * 3600 + blh1.L;

            B2 = blh2.B <= 90.0 * 3600 ? blh2.B : 90.0 * 3600;
            B2 = blh2.B >= -90.0 * 3600 ? blh2.B : -90.0 * 3600;
            L2 = blh2.L <= 180.0 * 3600 ? blh2.L : blh2.L - 360.0 * 3600;
            L2 = blh2.L >= -180.0 * 3600 ? blh2.L : 360.0 * 3600 + blh2.L;
            B1 = B1 / Angle.Rou;
            L1 = L1 / Angle.Rou;
            B2 = B2 / Angle.Rou;
            L2 = L2 / Angle.Rou;
            #endregion
            if (B1 == B2 && L1 == L2)
            {
                A1 = 0; A2 = 0; S = 0;
            }
            else
            {
                if (Math.Abs(B1) == 0.5 * Math.PI)
                    L1 = L2;
                if (Math.Abs(B2) == 0.5 * Math.PI)
                    L2 = L1;
                L = L2 - L1;
                while (L > Math.PI)
                    L -= 2.0 * Math.PI;
                while (L < -Math.PI)
                    L += 2.0 * Math.PI;
                //计算规划纬度u1,u2
                double u1 = 0;
                double u2 = 0;
                if (Math.Abs(B1) == 0.5 * Math.PI)
                    u1 = B1;
                else
                    u1 = Math.Atan((1 - blh1.RefEllipsoid.f) * Math.Tan(B1));
                if (Math.Abs(B2) == 0.5 * Math.PI)
                    u2 = B2;
                else
                    u2 = Math.Atan((1 - blh1.RefEllipsoid.f) * Math.Tan(B2));
                //迭代计算球面经差
                #region
                double w = 0;
                double delta = 0;
                double dw1 = 0;
                double dw2 = 0;
                double deltam = 0;
                double x = 0;
                double y = 0;
                double z = 0;
                double un = 0;
                do
                {
                    w = L + dw1;
                    while (w > Math.PI)
                        w -= 2.0 * Math.PI;
                    while (w <= -Math.PI)
                        w += 2.0 * Math.PI;
                    x = Math.Cos(w) * Math.Cos(u1) * Math.Cos(u2) + Math.Sin(u1) * Math.Sin(u2);
                    delta = Math.Acos(x);
                    //计算临界点的规划纬度un

                    if (Math.Abs(x) == 1)
                        y = 0;
                    else if (Math.Abs(x) != 1 && u1 == 0 && u2 == 0)
                        y = 1;
                    else
                        y = Math.Cos(u1) * Math.Abs(Math.Cos(u2) * Math.Sin(Math.Abs(w)) / Math.Sin(delta));
                    un = Math.Acos(y);
                    //计算大地线中心角距

                    if (y == 1)
                        z = x;
                    else
                        z = Math.Cos(delta) - 2.0 * Math.Sin(u1) * Math.Sin(u2) / (Math.Sin(un) * Math.Sin(un));
                    deltam = Math.Acos(z);
                    //计算经差改正项
                    double v = blh1.RefEllipsoid.f * Math.Sin(un) * Math.Sin(un) / 4.0;
                    double K3 = v * (1 + blh1.RefEllipsoid.f + blh1.RefEllipsoid.f * blh1.RefEllipsoid.f - v * (3 + 7 * blh1.RefEllipsoid.f - 13 * v));
                    double d = (1 - K3) * blh1.RefEllipsoid.f * Math.Cos(un) * (delta + K3 * Math.Sin(delta) * (Math.Cos(deltam) + K3 * Math.Cos(delta) * Math.Cos(2.0 * deltam)));
                    if (L == 0)
                        dw2 = 0;
                    else if (Math.Abs(L) == Math.PI)
                        dw2 = 0;
                    else if (L > 0)
                        dw2 = d;
                    else
                        dw2 = -d;
                    double temp = dw2;
                    dw2 = dw1;
                    dw1 = temp;
                } while (Math.Abs(dw1 - dw2) > 0.00000000000001);

                #endregion
                //计算大地线长s:
                double t = blh1.RefEllipsoid.seconde2 * Math.Sin(un) * Math.Sin(un) / 4.0;
                double K1 = 1 + t * (1 - t * (3 - t * (5 - 11 * t)) / 4.0);
                double K2 = t * (1 - t * (2 - t * (37 - 94 * t) / 8.0));
                double Ddelta = K2 * Math.Sin(delta) * (Math.Cos(deltam) + K2 * (Math.Cos(delta) * Math.Cos(2 * deltam) + K2 * (1 + 2 * Math.Cos(2 * delta)) * Math.Cos(3 * deltam) / 6) / 4);
                S = K1 * blh1.RefEllipsoid.b * (delta - Ddelta);
                //计算正反方位角
                if (L == 0 || Math.Abs(L) == Math.PI)
                {
                    if (L == 0 && B1 < B2)
                        A1 = 0;
                    if (L == 0 && B1 < B2)
                        A1 = Math.PI;
                    if (Math.Abs(L) == Math.PI && (Math.PI - B1 - B2) < (Math.PI + B1 + B2))
                        A1 = 0;
                    if (Math.Abs(L) == Math.PI && (Math.PI - B1 - B2) >= (Math.PI + B1 + B2))
                        A1 = Math.PI;
                    if (L == 0 && B1 > B2)
                        A2 = 0;
                    if (L == 0 && B1 <= B2)
                        A2 = Math.PI;
                    if (Math.Abs(L) == Math.PI)
                        A2 = A1;
                }
                else if (B1 == 0 && B2 == 0)
                {
                    if (L > 0)
                        A1 = 0.5 * Math.PI;
                    if (L <= 0)
                        A1 = 1.5 * Math.PI;
                    if (L < 0)
                        A2 = 0.5 * Math.PI;
                    if (L >= 0)
                        A2 = 1.5 * Math.PI;
                }
                else
                {
                    double m = Math.Cos(u2) * Math.Sin(w) / (Math.Cos(u1) * Math.Sin(u2) - Math.Sin(u1) * Math.Cos(u2) * Math.Cos(w));
                    double n = Math.Cos(u1) * Math.Sin(w) / (Math.Cos(u1) * Math.Sin(u2) * Math.Cos(w) - Math.Sin(u1) * Math.Cos(u2));
                    double sinw = Math.Sin(w);
                    if (sinw > 0 && m > 0)
                        A1 = Math.Atan(m);
                    if (sinw > 0 && m < 0)
                        A1 = Math.PI + Math.Atan(m);
                    if (sinw < 0 && m > 0)
                        A1 = Math.PI + Math.Atan(m);
                    if (sinw < 0 && m < 0)
                        A1 = 2.0 * Math.PI + Math.Atan(m);

                    if (sinw > 0 && n < 0)
                        A2 = 2.0 * Math.PI + Math.Atan(n);
                    if (sinw < 0 && n > 0)
                        A2 = 2.0 * Math.PI + Math.Atan(n);
                    if (sinw < 0 && n < 0)
                        A2 = Math.PI + Math.Atan(n);
                    if (sinw > 0 && n > 0)
                        A2 = Math.PI + Math.Atan(n);
                }
            }
            while (A1 >= 2.0 * Math.PI)
                A1 -= 2.0 * Math.PI;
            while (A1 < 0)
                A1 += 2.0 * Math.PI;
            while (A2 >= 2.0 * Math.PI)
                A2 -= 2.0 * Math.PI;
            while (A2 < 0)
                A2 += 2.0 * Math.PI;
            ang1 = new Angle(A1 * Angle.Rou);
            ang2 = new Angle(A2 * Angle.Rou);
            return S;
        }

        /// <summary>
        /// 计算某一纬度处子午圈曲率半径M
        /// </summary>
        /// <param name="refEll">椭球参数</param>
        /// <param name="b">纬度，以秒为单位</param>
        /// <returns>纬度为B的子午圈曲率半径</returns>
        public static double MeridianCircleRadius(ReferenceEllipsoid refEll, double b)
        {
            double res = 0;
            b = b / Angle.Rou;
            double sinB = Math.Sin(b);
            res = refEll.a * (1 - refEll.firste2) * Math.Pow(1 - refEll.firste2 * sinB * sinB, -1.5);
            return res;
        }
        /// <summary>
        /// 计算某一纬度处卯酉圈曲率半径N
        /// </summary>
        /// <param name="refEllipsoid">椭球参数</param>
        /// <param name="b">纬度，以秒为单位</param>
        /// <returns>纬度为B的卯酉圈曲率半径</returns>
        public static double UnitaryCircleRadius(ReferenceEllipsoid refEll, double b)
        {
            double res = 0;
            b = b / Angle.Rou;
            double sinB = Math.Sin(b);
            res = refEll.a * Math.Pow(1 - refEll.firste2 * sinB * sinB, -0.5);
            return res;
        }
        /// <summary>
        /// 计算某一纬度处平均曲率半径R
        /// </summary>
        /// <param name="refEll">椭球参数</param>
        /// <param name="b">纬度，以秒为单位</param>
        /// <returns>纬度为B的平均曲率半径</returns>
        public static double AverageCircleRadius(ReferenceEllipsoid refEll, double b)
        {
            double res = 0;
            double m = MeridianCircleRadius(refEll, b);
            double n = UnitaryCircleRadius(refEll, b);
            res = Math.Sqrt(m * n);
            return res;
        }
        /// <summary>
        /// 任意方向上的曲率半径
        /// </summary>
        /// <param name="refEll">参考椭球</param>
        /// <param name="B">纬度</param>
        /// <param name="A">方位角</param>
        /// <returns>该方向的曲率半径</returns>
        public static double AnyDiretionRadius(ReferenceEllipsoid refEll, Angle B, Angle A)
        {
            double N = UnitaryCircleRadius(refEll, B.Seconds);
            double cosB = Math.Cos(B.Radian);
            double cosA = Math.Cos(A.Radian);
            return N / (1 + refEll.seconde2 * cosB * cosB * cosA * cosA);
        }

        /// <summary>
        /// 计算从赤道开始到任意纬度B的平行圈之间的弧长
        /// </summary>
        /// <param name="refEll">椭球参数</param>
        /// <param name="B">纬度，以秒为单位</param>
        /// <returns>弧长，以米为单位</returns>
        public static double LengthOfMeridian(ReferenceEllipsoid refEll, double gB)
        {
            double res = 0;
            if (gB < 0)
                gB = -gB;
            gB = gB / Angle.Rou;
            double sinB = Math.Sin(gB);
            double cosB = Math.Cos(gB);
            double A1 = 1 + (3 / 4.0) * refEll.firste2 + (45 / 64.0) * Math.Pow(refEll.firste2, 2.0) + (175 / 256.0) * Math.Pow(refEll.firste2, 3.0) + (11025 / 16384.0) * Math.Pow(refEll.firste2, 4.0) + (43659 / 65536.0) * Math.Pow(refEll.firste2, 5.0) + (693693 / 1048576.0) * Math.Pow(refEll.firste2, 6.0);
            double B1 = (3 / 8.0) * refEll.firste2 + (15 / 32.0) * Math.Pow(refEll.firste2, 2.0) + (525 / 1024.0) * Math.Pow(refEll.firste2, 3.0) + (2205 / 4096.0) * Math.Pow(refEll.firste2, 4.0) + (72765 / 131072.0) * Math.Pow(refEll.firste2, 5.0) + (297297 / 524288.0) * Math.Pow(refEll.firste2, 6.0);
            double C1 = (15 / 256.0) * Math.Pow(refEll.firste2, 2.0) + (105 / 1024.0) * Math.Pow(refEll.firste2, 3.0) + (2205 / 16384.0) * Math.Pow(refEll.firste2, 4.0) + (10395 / 65536.0) * Math.Pow(refEll.firste2, 5.0) + (1486486 / 8388608.0) * Math.Pow(refEll.firste2, 6.0);
            double D1 = (35 / 3072.0) * Math.Pow(refEll.firste2, 3.0) + (105 / 4096.0) * Math.Pow(refEll.firste2, 4.0) + (10395 / 262144.0) * Math.Pow(refEll.firste2, 5.0) + (55055 / 1048576.0) * Math.Pow(refEll.firste2, 6.0);
            double E1 = (315 / 131072.0) * Math.Pow(refEll.firste2, 4.0) + (3465 / 524288.0) * Math.Pow(refEll.firste2, 5.0) + (99099 / 8388608.0) * Math.Pow(refEll.firste2, 6.0);
            double F1 = (693 / 1310720.0) * Math.Pow(refEll.firste2, 5.0) + (9009 / 524288.0) * Math.Pow(refEll.firste2, 6.0);
            double G1 = (1001 / 8388608.0) * Math.Pow(refEll.firste2, 6.0);
            res = refEll.a * (1 - refEll.firste2) * (A1 * gB - B1 * Math.Sin(2 * gB) + C1 * Math.Sin(4 * gB) - D1 * Math.Sin(6 * gB) + E1 * Math.Sin(8 * gB) - F1 * Math.Sin(10 * gB) + G1 * Math.Sin(12 * gB));
            return res;
        }
        /// <summary>
        /// 由子午线弧长计算该点的纬度
        /// </summary>
        /// <param name="refEll">椭球参数</param>
        /// <param name="gB">弧长，以米为单位</param>
        /// <returns>纬度，以秒为单位</returns>
        public static double LengthOfMeridianAntiCaculate(ReferenceEllipsoid refEll, double gx)
        {
            double res = 0;
            double b0 = gx * Math.Sqrt(1 - refEll.firste2) / refEll.a;
            double b1 = 0;
            do
            {
                double x1 = LengthOfMeridian(refEll, b0 * Angle.Rou);
                double deltax = gx - x1;
                double deltabeta = deltax / MeridianCircleRadius(refEll, b0 * Angle.Rou);
                double deltab = deltabeta - 1.5 * refEll.firste2 * ((1 + refEll.firste2 * Math.Sin(b0) * Math.Sin(b0)) * Math.Sin(2 * b0) * deltabeta * deltabeta / 2.0 + Math.Cos(2 * b0) * deltabeta * deltabeta * deltabeta / 3.0);
                b1 = b0 + deltab;
                double temp = b0;
                b0 = b1;
                b1 = temp;
            } while (Math.Abs(b0 - b1) > 0.000000000001);
            res = b0 * Angle.Rou;
            return res;
        }

        /// <summary>
        /// 由大地坐标计算子午线收敛角
        /// </summary>
        /// <param name="blh">大地坐标</param>
        /// <returns>子午线收敛角</returns>
        public static double MeridianConvergenceByBLH(GeodeticPoint blh)
        {
            double res = 0;
            double centerL = (blh.L + 3600 * 360.0) % (3600 * 360.0);
            int beltNo = (int)(centerL / (6.0 * 3600)) + 1;
            centerL = 6.0 * 3600 * beltNo - 3.0 * 3600;
            double l = blh.L - centerL;
            l = l / Angle.Rou;
            double b = blh.B / Angle.Rou;
            double sinB = Math.Sin(b);
            double cosB = Math.Cos(b);
            double t = Math.Tan(b);
            double yita2 = blh.RefEllipsoid.seconde2 * cosB * cosB;
            res = sinB * l + (1 / 3.0) * sinB * cosB * cosB * Math.Pow(l, 3) * (1 + 3 * yita2 + 2 * yita2 * yita2) + (1 / 15.0) * sinB * Math.Pow(cosB, 4) * Math.Pow(l, 5) * (2 - t * t);
            return res * Angle.Rou;
        }
        /// <summary>
        /// 由高斯平面坐标计算子午线收敛角
        /// </summary>
        /// <param name="xyh">高斯平面直角坐标</param>
        /// <returns>子午线收敛角</returns>
        public static double MeridianConvergenceByxy(GeodeticPoint xyh)
        {
            double bf = LengthOfMeridianAntiCaculate(xyh.RefEllipsoid, xyh.x);
            double y = xyh.y % 1000000 - 500000;
            double tf = Math.Tan(bf / Angle.Rou);
            double yitaf = Math.Sqrt(xyh.RefEllipsoid.seconde2) * Math.Cos(bf / Angle.Rou);
            double nf = UnitaryCircleRadius(xyh.RefEllipsoid, bf);
            return Angle.Rou * (y * tf / nf - Math.Pow(y, 3) * tf * (1 + tf * tf - yitaf * yitaf) / (3 * Math.Pow(nf, 3)) + Math.Pow(y, 5) * tf * (2 + 5 * tf * tf + 3 * Math.Pow(yitaf, 4)) / (15 * Math.Pow(nf, 5)));
        }
        /// <summary>
        /// 方向改化计算
        /// </summary>
        /// <param name="xyh1">高斯平面点1</param>
        /// <param name="xyh2">高斯平面点2</param>
        /// <returns>方向改化值</returns>
        public static double CurvatureCorrection(GeodeticPoint xyh1, GeodeticPoint xyh2)
        {
            double res = 0;
            double b1 = LengthOfMeridianAntiCaculate(xyh1.RefEllipsoid, xyh1.x);
            double b2 = LengthOfMeridianAntiCaculate(xyh2.RefEllipsoid, xyh2.x);
            double y1 = xyh1.y % 1000000 - 500000.0;
            double y2 = xyh2.y % 1000000 - 500000.0;
            double rm = AverageCircleRadius(xyh1.RefEllipsoid, (b1 + b2) / 2.0);
            res = -Angle.Rou * (xyh2.x - xyh1.x) * (2 * y1 + y2) / (6.0 * rm * rm);
            return res;
        }

        /// <summary>
        /// 由某方向的大地方位角、垂直角和垂线偏差计算垂线偏差改正
        /// </summary>
        /// <param name="kesi">垂线偏差子午线分量</param>
        /// <param name="yita">垂线偏差卯酉分量</param>
        /// <param name="ang">该方向上大地方位角</param>
        /// <param name="z">该方向的天顶距</param>
        /// <returns>垂线偏差改正 delta1</returns>
        public static double VerticalDeflectionCorrection(double kesi, double yita, Angle ang, Angle z)
        {
            double sita = Math.Atan2(yita, kesi) *Angle.Rou;
            sita = sita < 0 ? sita + 3600 * 360 : sita;
            if (kesi == 0 && yita == 0)
                return 0;
            if (z.Seconds == 90 * 3600.0)
                return 0;
            double sinA = Math.Sin(ang.Radian);
            double cosA = Math.Cos(ang.Radian);
            double tanalfa = Math.Tan(Math.PI / 2 - z.Radian);
            return (yita * cosA - kesi * sinA) * tanalfa;

        }
        /// <summary>
        /// 标高差改正计算delta2
        /// </summary>
        /// <param name="blh">照准点大地坐标</param>
        /// <param name="ang">照准点大地方位角</param>
        /// <returns>标高差改正数delta2,秒为单位</returns>
        public static double ElevationDifferenceCorrection(GeodeticPoint blh, Angle ang)
        {
            if (blh.H == 0)
                return 0;
            if ((ang.Seconds % 90.0 * 3600) == 0)
                return 0;
            double m = MeridianCircleRadius(blh.RefEllipsoid, blh.B);
            double cosB = Math.Cos(blh.B / Angle.Rou);
            double sin2A = Math.Sin(2 * ang.Radian);
            return blh.RefEllipsoid.firste2 * blh.H * Angle.Rou * cosB * cosB * sin2A / (2 * m);
        }
        /// <summary>
        /// 截面差改正delta3
        /// </summary>
        /// <param name="blh">测站点坐标</param>
        /// <param name="ang">大地方位角</param>
        /// <param name="s">法截线长</param>
        /// <returns>截面差改正delta3，秒为单位</returns>
        public static double NormalSectionCorrection(GeodeticPoint blh, Angle ang, double s)
        {
            if ((ang.Seconds % 90.0 * 3600) == 0)
                return 0;
            double n = UnitaryCircleRadius(blh.RefEllipsoid, blh.B);
            double cosB = Math.Cos(blh.B / Angle.Rou);
            double sin2A = Math.Sin(2 * ang.Radian);
            return blh.RefEllipsoid.firste2 * s * s * cosB * cosB * sin2A / (12 * n * n);
        }


        /// <summary>
        /// 由某方向的大地方位角、垂直角和垂线偏差计算观测天顶距改正
        /// </summary>
        /// <param name="kesi">垂线偏差子午线分量</param>
        /// <param name="yita">垂线偏差卯酉分量</param>
        /// <param name="ang">该方向上大地方位角</param>
        /// <param name="z">该方向的天顶距</param>
        /// <returns>归算后的天顶距，以秒为单位</returns>
        public static double ZenithCorrection(double kesi, double yita, Angle ang, Angle z)
        {
            double sinA = Math.Sin(ang.Radian);
            double cosA = Math.Cos(ang.Radian);
            return z.Seconds + kesi * cosA + yita * sinA;
        }
        /// <summary>
        /// 地面观测长度归算至椭球面
        /// </summary>
        /// <param name="refEll">椭球参数</param>
        /// <param name="h1">测站点大地高</param>
        /// <param name="h2">照准点大地高</param>
        /// <param name="d">已知斜距</param>
        /// <param name="b">测站点大地纬度</param>
        /// <param name="a">测距边大地方位角</param>
        /// <returns>椭球面大地线长</returns>
        public static double EllipsoidDistanceCorrection(ReferenceEllipsoid refEll, double h1, double h2, double d, Angle b, Angle a)
        {
            double ra = AnyDiretionRadius(refEll, b, a);
            double hm = (h1 + h2) / 2.0;
            double dd = Math.Sqrt(d * d - (h2 - h1) * (h2 * h1));
            double sin2B = Math.Sin(2 * b.Radian);
            double cosA = Math.Cos(a.Radian);
            return dd * ra / (ra + hm) + d * d * d / (24 * ra * ra) + 1.25E-16 * hm * d * d * sin2B * cosA;
        }
        /// <summary>
        /// 由近似坐标计算高斯距离改正值
        /// </summary>
        /// <param name="xyh1">起点近似坐标</param>
        /// <param name="xyh2">终点近似坐标</param>
        /// <param name="s">大地线长</param>
        /// <returns>大地线化为高斯平面距离的改正值</returns>
        public static double GaussDistanceCorrection(GeodeticPoint xyh1, GeodeticPoint xyh2, double s)
        {
            double bm = LengthOfMeridianAntiCaculate(xyh1.RefEllipsoid, (xyh1.x + xyh2.x) / 2.0);
            double rm = AverageCircleRadius(xyh1.RefEllipsoid, bm);
            double ym = (xyh1.y % 1000000 + xyh2.y % 1000000 - 1000000) / 2;
            double dy = xyh2.y - xyh1.y;
            return s * (0.5 * ym * ym / (rm * rm) + dy * dy / (24 * rm * rm) + 1 / 24.0 * (Math.Pow(ym / rm, 4)));
        }
    }
}
