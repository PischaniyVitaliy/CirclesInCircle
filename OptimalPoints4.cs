using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Cureos.Numerics;
/* Date: 18-02-2017
 * Круги в круге 
 * countCircle ~ можно задавать количество степеней в ограничении x[i]^k = r[i]^k 
 * Исходные данные наборами кругов одинакового радиуса - l
 ******************************** */
namespace hs071_cs
{
  // Классы задач 
  public enum Restriction
  {
    Power /*степенные ограничения*/,
    Lineal /*линейные ограничения*/,
    AlphaPower
  }

  public class OptimalPoints4
  {
    public readonly int _n;
    private readonly Restriction _task;
    public int _m;
    public int _nele_jac;
    public int _nele_hess;
    public double[] _x_L;
    public double[] _x_U;
    public double[] _g_L;
    public double[] _g_U;

    public double[] X { get; private set; } // Для решателя массив х, в него входят все координаты и радиусы. Т.е. переменные

    public readonly int count; // количество кругов //альтернативная задача => this.height = height;
    public double[] radius;
    private int[] _groupNumber; // указывает номер группы в котором находится круг
    private int _limitCount; // количество степеней в ограничении x[i]^k = r[i]^k
    private readonly int _groupCount;
    private readonly int[] _elemGroupCount;

    // Начальная точка - вычисляется программно на основе радиусов. Метод: FirstValue(r).
    // countCircle - степень в ограничениях 
    public OptimalPoints4(double[] r, int[] groupNum, Restriction task = Restriction.Power, int countCircle = 1, double[] x = null, double[] y = null)
    {
      _task = task; // передаём тип ограничений решаемой в задаче
      count = r.Length; // задаём количетво кругов
      _n = count * 3 + 1; // количество переменных в векторе 
      radius = new double[count];
      int it = 0;
      foreach (var item in r)
        radius[it++] = item;
      X = new double[_n];
      // Если передаём исходные x и y, то инициализируем ими, иначе генерим
      if (x != null && y != null)
      {
        if (x.Length != y.Length || x.Length != r.Length || x.Length != groupNum.Length) throw new Exception();
        for (int i = 0; i < count; ++i)
        {
          X[2 * i] = x[i];
          X[2 * i + 1] = y[i];
          X[2 * count + i] = r[i];
        }
      }
      else
      {
        if (r.Length != groupNum.Length) throw new Exception();
        // генеририруем R из массивов Х и У
        x = new double[count];
        y = new double[count];
        double l = 0;
        double d = 0.2;
        for (int i = 0; i < count; ++i)
        {
          
          X[2 * i] = x[i] = l + r[i];
          l += 2 * r[i];
          X[2 * i + 1] = y[i] = d *= -1;
          X[2 * count + i] = r[i];
          /*
          X[2 * i] = 0;
          X[2 * i + 1] = 0;
          X[2 * count + i] = 0;
          */
        }
        X[_n - 1] = l;
      }
      switch (_task)
      {
        case Restriction.Lineal:
          break;
        case Restriction.Power:
          _limitCount = countCircle;
          break;
      }
      /*    Вычисляем количество групп
      * *************************************************************************************/
      _groupCount = 0;
      _groupNumber = new int[count];
      for (int i = 0; i < count; ++i)
      {
        _groupNumber[i] = groupNum[i];
        if (_groupCount < groupNum[i]) _groupCount = groupNum[i];
      }
      // находим сколько в каждой группе эл-ов, т.к. они не подряд 
      // ещё нулевая группа (в которой фиксированные радиусы)
      _groupCount++;
      _elemGroupCount = new int[_groupCount];
      for (int i = 0; i < count; ++i)
      {
        _elemGroupCount[groupNum[i]]++;
      }

      /*    Ограничения переменных
      * *************************************************************************************/
      _x_L = new double[_n];
      double max = radius[radius.Length - 1];

      for (var i = 0; i < radius.Length; ++i)
      {
        _x_L[2 * i] = _x_L[2 * i + 1] = -DiametrSum(radius) + radius[i];
        _x_L[2 * count + i] = 0;
        //if (max < radius[i]) max = radius[i];
      }
      _x_L[_n - 1] = max;

      _x_U = new double[_n];
      for (var i = 0; i < count; ++i)
      {
        _x_U[2 * i] = _x_U[2 * i + 1] = DiametrSum(radius) - radius[i];
        _x_U[2 * count + i] = radius[radius.Length - 1]; // max Radius
      }
      _x_U[_n - 1] = DiametrSum(radius) / 2; //Ipopt.PositiveInfinity;

      /*    Огрничения
       **************************************************************************************/
      // не нулевых ограничений в Якобиане по первым двум группам ограничений
      _nele_jac = 0;
      _m = 0;
      // (R-r[i])^2-x[i]^2-y[i]^2 >= 0
      _nele_jac += 4 * count;
      _m += count;
      // (x[i]-x[j])^2 + (y[i]-y[j])^2 - (r[i]-r[j])^2 >=0
      _nele_jac += 6 * count * (count - 1) / 2;
      _m += count * (count - 1) / 2;

      switch (_task)
      {
        // ~~~~~~~~~~~~~~~~~~~~~~~~ Tasks.Lineal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case Restriction.Lineal:
          var countLineal = (int)Math.Pow(2, count);
          _m += countLineal - 1; // все линейные ограничения
          _m += 1;               // sum{r[i]^2} = sum{vr[i]^2} 

          _g_L = new double[_m];
          _g_U = new double[_m];

          for (var i = 0; i < _m - 2; i++)
          {
            _g_L[i] = 0;
            _g_U[i] = Ipopt.PositiveInfinity;
          }
          // На равенство
          for (var i = _m - 2; i < _m; i++)
          {
            _g_L[i] = 0;
            _g_U[i] = 0;
          }

          int oneCount = 0; // количество едениц
          for (int i = 1; i < countLineal; ++i)
          {
            var s = Convert.ToString(i, 2);
            for (int k = 0; k < s.Length; ++k)
              if (s[k] == '1')
                oneCount++;
          }

          // всевозможные комбинации линейных ограничений
          _nele_jac += oneCount;
          // сумма квадратов sum{r[i]^2} = sum{vr[i]^2} 
          _nele_jac += count;

          _nele_hess = 0;
          break;
        // ~~~~~~~~~~~~~~~~~~~~~~~~ Tasks.Power ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case Restriction.Power:
          _m += count;
          _g_L = new double[_m];
          _g_U = new double[_m];

          // Для любых ограничений на неравенство
          for (var i = 0; i < _m - count; i++)
          {
            _g_L[i] = 0;
            _g_U[i] = Ipopt.PositiveInfinity;
          }

          // Дополнительные ограничения sum{x[i]} - sum{vx[i]} = 0
          for (var i = _m - count; i < _m; i++)
          {
            _g_L[i] = 0;
            _g_U[i] = 0;
          }

          // не нулевых ограничений в Якобиане (3-я группа ограничений)
          int _nele_jac_g3 = 0;
          for (var i = 0; i < _groupCount; ++i)
            _nele_jac_g3 += _elemGroupCount[i] * _elemGroupCount[i];
          _nele_jac += _nele_jac_g3;

          break;
      } // switch (_task)
    } // End_Конструктор 

    public double DiametrSum(double[] radius)
    {
      var sum = 0.0;
      foreach (var rad in radius)
      {
        sum += 2 * rad;
      }
      return sum;
    }

    public double[] FirstValue(double[] r)
    {
      double[] x = new double[count * 3 + 1];
      double l = 0;
      double d = 0.2;
      for (int i = 0; i < count; ++i)
      {
        x[2 * i] = l + r[i];
        l += 2 * r[i];
        x[2 * i + 1] = d *= -1;
        //x[2 * count + i] = r[i];
        x[2 * count + i] = 20;
      }
      x[_n - 1] = l;
      return x;
    }

    // Инициализация значений всего массива _х
    public double[] InitialValuesManually(double[] InitialValues)
    {
      double[] x = new double[count * 3 + 1];
      for (int i = 0; i < count; ++i)
      {
        x[2 * i] = InitialValues[2 * i];
        x[2 * i + 1] = InitialValues[2 * i + 1];
        x[2 * count + i] = InitialValues[2 * count + i];
      }
      x[_n - 1] = InitialValues[_n - 1];
      return x;
    }
    public bool eval_f(int n, double[] x, bool new_x, out double obj_value)
    {
      obj_value = x[_n - 1];

      return true;
    }

    public bool eval_grad_f(int n, double[] x, bool new_x, double[] grad_f)
    {
      for (var i = 0; i < _n - 1; i++)
      {
        grad_f[i] = 0;
      }
      grad_f[_n - 1] = 1;
      return true;
    }

    string DecToBase(int num_value, int base_value, int max_bit = 64)
    {
      var dec_base = 10;
      var hexchars = new[] { 'A', 'B', 'C', 'D', 'E', 'F' };
      var result = string.Empty;
      var result_array = new int[max_bit];

      for (/* nothing */; num_value > 0; num_value /= base_value)
      {
        int i = num_value % base_value;
        result_array[--max_bit] = i;
      }

      for (int i = 0; i < result_array.Length; i++)
      {
        if (result_array[i] >= dec_base)
        {
          result += hexchars[(int)result_array[i] % dec_base].ToString();
        }
        else
        {
          result += result_array[i].ToString();
        }
      }

      //result = result.TrimStart(new char[] {'0'});
      return result;
    }

    public bool eval_g(int n, double[] x, bool new_x, int m, double[] g)
    {
      int kk = 0;
      // from 0 to count-1
      for (var i = 0; i < count; ++i)
        g[kk++] = Math.Pow((x[_n - 1] - x[count * 2 + i]), 2) - x[2 * i] * x[2 * i] - x[2 * i + 1] * x[2 * i + 1];
      //g[i] = x[_n - 1] - radius[i] - Math.Sqrt(x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1]);
      // kk = count
      // from count to count*(count-1)/2 - 1
      for (var i = 0; i < count - 1; ++i)
      {
        for (var j = i + 1; j < count; ++j)
          g[kk++] = Math.Pow((x[2 * i] - x[2 * j]), 2) + Math.Pow((x[2 * i + 1] - x[2 * j + 1]), 2) - Math.Pow((x[count * 2 + i] + x[count * 2 + j]), 2);
      }
      // kk=count+count*(count-1)/2 = count + count^2/2 - count/2 = (count + count^2)/2
      switch (_task)
      {
        // ~~~~~~~~~~~~~~~~~~~~~~~~ Tasks.Lineal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case Restriction.Lineal:
          var sumRadius = new double[count];
          for (int i = 0; i < count; ++i)
            for (int j = 0; j <= i; ++j)
              sumRadius[i] += radius[j];
          var countLineal = Math.Pow(2, count);
          var oneCount = 0;
          for (int a = 1; a < countLineal; ++a)
          {
            var bitStr = DecToBase(a, 2, count);
            oneCount = 0;
            for (int b = 0; b < bitStr.Length; ++b)
            {
              if (bitStr[b] == '1')
              {
                g[kk] += x[2 * count + b];
                oneCount++;
              }
            }
            g[kk++] -= sumRadius[oneCount - 1]; // количество едениц, но индексы с 0
          }
          // сумма квадратов sum{r[i]^2} = sum{vr[i]^2} 
          g[kk] = 0;
          for (int i = 0; i < count; ++i)
          {
            g[kk] += Math.Pow(x[2 * count + i],2);
            g[kk] -= Math.Pow(radius[i], 2);
          }
          kk++;

          break;
        // ~~~~~~~~~~~~~~~~~~~~~~~~ Tasks.Power ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case Restriction.Power:
          double sumX = 0, sumVX = 0;
          int rank = 1;
          for (int k = 1; k < _groupCount; ++k) // до количества групп
          {
            rank = 1; // степень уравнений ограничения
            sumX = 0;
            sumVX = 0;
            for (int a = 0; a < _elemGroupCount[k]; ++a) // до количества эл-ов в каждой группе
            {
              for (int i = 0; i < count; ++i)
              {
                if (k == _groupNumber[i])
                {
                  sumX += Math.Pow(x[2 * count + i], rank);
                  sumVX += Math.Pow(radius[i], rank);
                }
              }
              g[kk++] = sumX - sumVX;
              rank++;
            }
          }
          break;
      }

      return true;
    }

    public bool eval_jac_g(int n, double[] x, bool new_x, int m, int nele_jac, int[] iRow, int[] jCol, double[] values)
    {
      if (values == null)
      {
        int kk = 0;
        var g = 0;
        // позиции по Х и У
        for (g = 0; g < count; ++g)
        {
          iRow[kk] = g;
          jCol[kk++] = _n - 1;
          iRow[kk] = g;
          jCol[kk++] = 2 * g;
          iRow[kk] = g;
          jCol[kk++] = 2 * g + 1;
          iRow[kk] = g;
          jCol[kk++] = count * 2 + g;
        }
        for (var i = 0; i < count - 1; ++i)
        {
          for (var j = i + 1; j < count; ++j)
          {
            iRow[kk] = g;
            jCol[kk++] = 2 * i;
            iRow[kk] = g;
            jCol[kk++] = 2 * j;
            iRow[kk] = g;
            jCol[kk++] = 2 * i + 1;
            iRow[kk] = g;
            jCol[kk++] = 2 * j + 1;
            iRow[kk] = g;
            jCol[kk++] = 2 * count + i;
            iRow[kk] = g++;
            jCol[kk++] = 2 * count + j;
          }
        }
        // g[3]
        switch (_task)
        {
          // ~~~~~~~~~~~~~~~~~~~~~~~~ Tasks.Lineal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          case Restriction.Lineal:
            var countLineal = Math.Pow(2, count); // кол-во ограничений
            for (int a = 1; a < countLineal; ++a)
            {
              var bitStr = DecToBase(a, 2, count);
              for (int b = 0; b < bitStr.Length; ++b)
              {
                /*
                if (bitStr[b] == '1')
                {
                  iRow[kk] = g;
                  jCol[kk++] = b * 2; // x[i]
                }
                */
                if (bitStr[b] == '1')
                {
                  iRow[kk] = g;
                  jCol[kk++] = 2 * count + b; // r[i]
                }
              }
              g++;
            }
            // sum{r[i]^2} = sum{vr[i]^2} 
            for (int a = 0; a < count; ++a)
            {
              iRow[kk] = g;
              jCol[kk++] = 2 * count + a;
            }
            break;
          // ~~~~~~~~~~~~~~~~~~~~~~~~ Tasks.Power ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          case Restriction.Power:
            for (var gr = 1; gr < _groupCount; ++gr) // перебираем группы
            {
              for (var line = 0; line < _elemGroupCount[gr]; ++line) // перебираем по строкам
              {
                for (var i = 0; i < count; ++i) // перебираем по столбцам
                {
                  if (gr == _groupNumber[i])
                  {
                    iRow[kk] = g;
                    jCol[kk++] = 2 * count + i;
                  }
                }
                g++;
              }
            }
            break;
        }

      }
      else
      {
        int kk = 0;
        for (int i = 0; i < count; ++i)
        {
          values[kk++] = 2 * (x[_n - 1] - x[2 * count + i]); // R0'
          values[kk++] = -2 * x[2 * i]; //X'
          values[kk++] = -2 * x[2 * i + 1]; //Y' 
          values[kk++] = -2 * (x[_n - 1] - x[2 * count + i]); //r'
          //values[kk++] = 1.0;
          //values[kk++] = -x[2 * i] / Math.Sqrt(x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1]);
          //values[kk++] = -x[2 * i + 1] / Math.Sqrt(x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1]);
        }
        for (var i = 0; i < count - 1; ++i)
        {
          for (var j = i + 1; j < count; ++j)
          {
            values[kk++] = 2 * (x[2 * i] - x[2 * j]); //X[i]'
            values[kk++] = -2 * (x[2 * i] - x[2 * j]); //X[j]'
            values[kk++] = 2 * (x[2 * i + 1] - x[2 * j + 1]); //Y[i]'
            values[kk++] = -2 * (x[2 * i + 1] - x[2 * j + 1]); //Y[j]'
            values[kk++] = -2 * (x[2 * count + i] + x[2 * count + j]); //R[i]'
            values[kk++] = -2 * (x[2 * count + i] + x[2 * count + j]); //R[j]'
          }
        }

        // g[3]
        switch (_task)
        {
          // ~~~~~~~~~~~~~~~~~~~~~~~~ Tasks.Lineal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          case Restriction.Lineal:
            var countLineal = Math.Pow(2, count); // кол-во ограничений
            for (int a = 1; a < countLineal; ++a)
            {
              var bitStr = DecToBase(a, 2, count);
              for (int b = 0; b < bitStr.Length; ++b)
              {
                if (bitStr[b] == '1')
                {
                  //values[kk++] = 1; // x[i]
                  values[kk++] = 1; //r[i]
                }
              }
            }
            // sum{r[i]^2} = sum{vr[i]^2} 
            for (int a = 0; a < count; ++a)
            {
              values[kk++] = 2 * x[2 * count + a]; //r[i]
            }
            break;
          // ~~~~~~~~~~~~~~~~~~~~~~~~ Tasks.Power ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          case Restriction.Power:
            int rank; // степень в уравненях
            for (var gr = 1; gr < _groupCount; ++gr) // перебираем группы
            {
              rank = 0;
              for (var line = 0; line < _elemGroupCount[gr]; ++line) // перебираем по строкам
              {
                rank++;
                for (var i = 0; i < count; ++i) // перебираем по столбцам
                {
                  if (gr == _groupNumber[i])
                  {
                    values[kk++] = rank * Math.Pow(x[2 * count + i], rank - 1);
                  }
                }
              }
            }
            break;
        } // *** switch (_task) ***
      }

      return true;
    }

    public bool eval_h(int n, double[] x, bool new_x, double obj_factor,
               int m, double[] lambda, bool new_lambda,
               int nele_hess, int[] iRow, int[] jCol,
               double[] values)
    {
      return false;
    }

    public override string ToString()
    {
      return "OptimalPoints4";
    }

  }
}

