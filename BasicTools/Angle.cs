using System;
using System.Collections.Generic;
using System.Text;

namespace TailaoZLY.BasicTools
{
    /// <summary>
    /// 格式化输入或输出一个角度
    /// 四种格式：度、弧度，秒、度分秒连写
    /// 接受四种格式输入，支持文本、浮点数输入
    /// 可以进行四种格式字符串输出
    /// </summary>
    public class Angle
    {
        const double Epsilon = 1.0E-12;
        public static double Rou = 206264.80624709635515647335733078;
        private double _seconds = 0;
        public double Seconds { get { return _seconds % (3600 * 1800); } }
        public double Radian { get { return _seconds / Rou; } }

        public Angle(string str, AngleStyle aStyle)
        {
            double dms = 0;
            if (double.TryParse(str, out dms))
                _seconds = new Angle(dms, aStyle).Seconds;
        }
        public Angle(String str)
        {
            _seconds = new Angle(str, AngleStyle.Dms).Seconds;
        }

        public Angle(double dms, AngleStyle aStyle)
        {
            if (aStyle == AngleStyle.Dms)
            {
                int deg = (int)(dms + Epsilon) / 10000;
                int min = (int)(dms - deg * 10000 + Epsilon) / 100;
                double sec = dms - deg * 10000 - min * 100;
                _seconds = deg * 3600 + min * 60 + sec;

            }
            else if (aStyle == AngleStyle.Deg)
            {
                _seconds = dms * 3600;
            }
            else if (aStyle == AngleStyle.Rad)
            {
                _seconds = dms * Rou;
            }
            else
            {
                _seconds = dms;
            }
        }
        public Angle(double dms)
        {
            _seconds = new Angle(dms, AngleStyle.Sec).Seconds;
        }

        public string ToString(AngleStyle astyle, int digits)
        {
            string res = "";
            double dms = (_seconds < 0) ? -_seconds : _seconds;
            string sign = (_seconds < 0) ? "-" : "";
            double dec = 0.0;
            if (astyle == AngleStyle.Dms)
            {
                int deg = (int)(dms / 3600 + Epsilon);
                int min = (int)(dms - deg * 3600 + Epsilon) / 60;
                int sec = (int)(dms - deg * 3600 - min * 60 + Epsilon);
                dec = dms - deg * 3600 - min * 60 - sec;
                res = sign + string.Format("{0}{1:00}{2:00}", deg, min, sec);
            }
            if (astyle == AngleStyle.Deg)
            {
                int deg = (int)(dms / 3600 + Epsilon);
                dec = dms / 3600.0 - deg;
                dec = Math.Round(dec, digits);
                res = sign + string.Format("{0}", deg);
            }
            if (astyle == AngleStyle.Rad)
            {
                int deg = (int)(dms / Rou + Epsilon);
                dec = dms / Rou - deg;
                res = sign + string.Format("{0}", deg);
            }
            System.Globalization.NumberFormatInfo provider = new System.Globalization.NumberFormatInfo
            {
                NumberDecimalDigits = digits
            };
            string str = dec.ToString("N", provider);
            str = str.Remove(0, 1);
            res += str;
            return res;
        }
        public string ToString(AngleStyle astyle)
        {
            return ToString(astyle, 8);
        }
        public override string ToString()
        {
            return ToString(AngleStyle.Dms, 8);
        }
        public static Angle operator +(Angle ang1, Angle ang2)
        {
            double value = ang1.Seconds + ang2.Seconds;
            return new Angle(value);
        }
        public static Angle operator -(Angle ang1, Angle ang2)
        {
            double value = ang1.Seconds - ang2.Seconds;
            return new Angle(value);
        }
    }
    public enum AngleStyle
    {
        Deg,
        Dms,
        Rad,
        Sec
    }
}
