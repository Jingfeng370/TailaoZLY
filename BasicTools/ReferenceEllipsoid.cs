using System;
using System.Collections.Generic;
using System.Text;

namespace TailaoZLY.BasicTools
{
    /// <summary>
    /// 参考椭球，有长半轴和扁率的倒数定义
    /// </summary>
    public class ReferenceEllipsoid
    {
        #region

        private double _a = 0;
        /// <summary>
        /// 长半轴
        /// </summary>
        public double a { get { return _a; }  }

        
        private double _alfa = 0;
        /// <summary>
        /// 扁率的倒数
        /// </summary>
        public double alfa { get { return _alfa; }  }

        
        private double _b = 0;
        /// <summary>
        /// 短半轴
        /// </summary>
        public double b { get { return _b; } }

        
        private double _c = 0;
        /// <summary>
        /// 极曲率半径
        /// </summary>
        public double c { get { return _c; } }

        
        private double _firste2 = 0;
        /// <summary>
        /// 第一偏心率的平方
        /// </summary>
        public double firste2 { get { return _firste2; } }

        
        private double _seconde2 = 0;
        /// <summary>
        /// 第二偏心率平方
        /// </summary>
        public double seconde2 { get { return _seconde2; } }

        
        private double _f = 0;
        /// <summary>
        /// 扁率
        /// </summary>
        public double f { get { return _f; } }

        public ReferenceEllipsoid(double a, double alfa)
        {
            _a = a;
            _b = a * (1 - alfa);
            _c = _a * _a / _b;
            _alfa = alfa;
            _firste2 = 1 - (_b * _b) / (_a * _a);
            _seconde2 = (a * a) / (b * b) - 1;
            _f = 1 - _b / _a;
        }
        #endregion//参考椭球相关参数
        public static readonly ReferenceEllipsoid WGS84 = new ReferenceEllipsoid(6378137, 1 / 298.257223563);
        public static readonly ReferenceEllipsoid CGCS2000 = new ReferenceEllipsoid(6378137, 1 / 298.257222101);
        public static readonly ReferenceEllipsoid Krassovsky = new ReferenceEllipsoid(6378245, 1 / 298.3);

        /// <summary>
        /// 通过字符串返回相应的椭球
        /// </summary>
        /// <param name="str">默认返回WGS84椭球，C代表CGCS2000椭球</param>
        /// <returns></returns>
        public static ReferenceEllipsoid GetRefEll(string str)
        {
            str.ToUpper();
            if (str == "C" || str == "CGCS2000")
                return CGCS2000;
            else if (str == "B" || str == "BJ54")
                return Krassovsky;
            else
                return WGS84;
        }
    }
}
