// Сообщения об ошибках ipopt: https://www.coin-or.org/Ipopt/documentation/node36.html
using System;
using Cureos.Numerics;
using System.IO;
using System.Diagnostics;

namespace hs071_cs
{
  public class Program
  {
    private static CircleTaskTypes taskType;
    private static Restriction restrictionType;
    private static Random rnd = new Random();
    public static void Main(string[] args)
    {
      taskType = CircleTaskTypes.DynamicRadius;
      restrictionType = Restriction.Power;
      if (taskType == CircleTaskTypes.FixedRadius) restrictionType = Restriction.FixedRadius;
      double[] r;
      int[] gr;
      double[] xNach;
      double[] yNach;
      OptimalPoints4 p;

      int colGenerateRadius = 10, maxRadius = 9;
      r = new double[colGenerateRadius];
      for (int i = 0; i < colGenerateRadius; ++i)
        r[i] = rnd.Next(1+maxRadius); // 1..maxRadius+1

      if (args.Length != 0)
      {
        OpenFromFile(args[0], out r, out gr, out restrictionType, out xNach, out yNach);
        p = new OptimalPoints4(r, gr, restrictionType);
        //p = new OptimalPoints4(r, gr, restrictionType, xNach, yNach);
      }
      else
      {
        r = new double[] { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1 };
        p = new OptimalPoints4(r, new int[] { 1, 1, 1, 1, 2, 2, 2 }, Restriction.Lineal);
      }
      // p = new OptimalPoints4(r, new int[] { 1, 1, 1, 1, 1, 2, 2, 2, 2, 2 }, restrictionType);
      // p = new OptimalPoints4(r, new int[] { 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 }, restrictionType);
      p = new OptimalPoints4(r, new int[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, restrictionType);


      var x = p.X;
      switch (taskType)
      {
        case CircleTaskTypes.FixedRadius:
          ShowData(p.X, r);
          break;
        case CircleTaskTypes.DynamicRadius:
          ShowData(p.X);
          break;
      }
      Console.WriteLine("\n\n<<<Нажмите любую клавишу>>>");
      Console.ReadKey();

      Stopwatch stopWatch = new Stopwatch();
      stopWatch.Start();
      /* allocate space for the initial point and set the values */
      IpoptReturnCode status;
      using (Ipopt problem = new Ipopt(p._n, p._x_L, p._x_U, p._m, p._g_L, p._g_U, p._nele_jac, p._nele_hess, p.eval_f, p.eval_g, p.eval_grad_f, p.eval_jac_g, p.eval_h))
      {
        /* Set some options.  The following ones are only examples,
           they might not be suitable for your problem. */
        // https://www.coin-or.org/Ipopt/documentation/node41.html#opt:print_options_documentation
        problem.AddOption("tol", 1e-7);
        problem.AddOption("mu_strategy", "adaptive");
        problem.AddOption("hessian_approximation", "limited-memory");
        problem.AddOption("output_file", p.ToString() + "_" + DateTime.Now.ToShortDateString() + "_"
          + DateTime.Now.Hour + "_" + DateTime.Now.Minute + ".txt");
        problem.AddOption("print_frequency_iter", 20);
        //problem.AddOption("file_print_level", 7); // 0..12
        problem.AddOption("file_print_level", 0);
        problem.AddOption("max_iter", 100000);
        problem.AddOption("print_level", 5); // 0<= value <= 12, default value is 5
#if INTERMEDIATE
                problem.SetIntermediateCallback(p.intermediate);
#endif
        /* solve the problem */
        double obj;
        status = problem.SolveProblem(x, out obj, null, null, null, null);
      }
      Console.WriteLine("{0}{0}Optimization return status: {1}{0}{0}", Environment.NewLine, status);
      stopWatch.Stop();

      switch (taskType)
      {
        case CircleTaskTypes.FixedRadius:
          ShowData(x, r);
		  //SaveToFile("result.txt", r, x, n);
          break;
        case CircleTaskTypes.DynamicRadius:
          ShowData(x);
		  
          break;
      }
      SaveToFile("Result.txt", x);
      Console.WriteLine("RunTime: " + getElapsedTime(stopWatch));
      Console.WriteLine("{0}Press <RETURN> to exit...", Environment.NewLine);
      Console.ReadLine();
    }

    /// <summary>
    /// Форматирует результат конвертирования времени запуска программы 
    /// </summary>
    /// <param name="Watch">объект Stopwatch</param>
    /// <returns>Строка- время чч:мм:сс.мс</returns>
    static string getElapsedTime(Stopwatch Watch)
    {
      // Get the elapsed time as a TimeSpan value.
      TimeSpan ts = Watch.Elapsed;
      // Format and display the TimeSpan value.
      return String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
    }

    /// <summary>
    /// Записывает результаты расчётов (для динамических радиусов)
    /// </summary>
    /// <param name="fileName">Путь к файлу</param>
    /// <param name="r">Радиусы всех кругов</param>
    /// <param name="x">Вектор оптимизируемый</param>
    /// <param name="sizeC">Количество кругов</param>
    static void SaveToFile(string fileName, double[] x, Restriction rType = Restriction.FixedRadius, double[] group = null)
    {
      StreamWriter sw = new StreamWriter(fileName, false);
      sw.WriteLine("# Тип ограничений");
      sw.WriteLine("# 0 = При фиксированных радиуасах нет дополнительных ограничений");
      sw.WriteLine("# 1 = Степенные");
      sw.WriteLine("# 2 = Линейные (многогранник)");
      int countCircle = (x.Length - 1) / 3;
      switch (restrictionType)
      {
        case Restriction.Lineal:
          sw.WriteLine("2");
          break;
        case Restriction.Power:
          sw.WriteLine("1");
          break;
      }
      sw.WriteLine("# Количество кругов");
      sw.WriteLine(countCircle);
      sw.WriteLine("# Радиусы");
      for (int i = 0; i < countCircle; i++)
      {
        sw.Write("{0} ", x[2 * countCircle + i]);
      }
      sw.WriteLine();
      sw.WriteLine("# Группы");
      if (group != null)
      {
        for (int i = 0; i < countCircle; i++)
        {
          sw.Write("{0} ", group[i]);
        }
      }
      else
      {
        for (int i = 0; i < countCircle; i++)
        {
          sw.Write("1 ");
        }
      }
      sw.WriteLine();
      sw.WriteLine("# Будут ли задаваться координаты х и у");
      sw.WriteLine("1");
      sw.WriteLine("# Координаты Х");
      for (int i = 0; i < countCircle; i++)
      {
        sw.Write("{0} ", x[2 * i]);
      }
      sw.WriteLine();
      sw.WriteLine("# Координаты У");
      for (int i = 0; i < countCircle; i++)
      {
        sw.Write("{0} ", x[2 * i + 1]);
      }
      sw.WriteLine();
      sw.WriteLine("# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ");
      sw.WriteLine("# Радиуc внешнего круга");
      sw.WriteLine("{0} ", x[3 * countCircle]);
      sw.Close();
    }

    /// <summary>
    /// Записывает результаты расчётов (для фиксированных радиусов)
    /// </summary>
    /// <param name="fileName">Путь к файлу</param>
    /// <param name="r">Радиусы всех кругов</param>
    /// <param name="x">Вектор оптимизируемый</param>
    /// <param name="sizeC">Количество кругов</param>
    static void SaveToFile(string fileName, double[] r, double[] x, Restriction rType = Restriction.FixedRadius, double[] group = null)
    {
      StreamWriter sw = new StreamWriter(fileName, false);
      sw.WriteLine("# Тип ограничений");
      sw.WriteLine("# 0 = При фиксированных радиуасах нет дополнительных ограничений");
      sw.WriteLine("# 1 = Степенные");
      sw.WriteLine("# 2 = Линейные (многогранник)");
      int countCircle = (x.Length - 1) / 3;
      switch (restrictionType)
      {
        case Restriction.Lineal:
          sw.WriteLine("2");
          break;
        case Restriction.Power:
          sw.WriteLine("1");
          break;
        case Restriction.FixedRadius:
          sw.WriteLine("0");
          break;
      }
      sw.WriteLine("# Количество кругов");
      sw.WriteLine(countCircle);
      sw.WriteLine("# Радиусы");
      for (int i = 0; i < countCircle; i++)
      {
        sw.Write("{0} ", r[i]);
      }
      sw.WriteLine();
      sw.WriteLine("# Группы");
      if (group != null)
      {
        for (int i = 0; i < countCircle; i++)
        {
          sw.Write("{0} ", group[i]);
        }
      }
      else
      {
        for (int i = 0; i < countCircle; i++)
        {
          sw.Write("1 ");
        }
      }
      sw.WriteLine();
      sw.WriteLine("# Будут ли задаваться координаты х и у");
      sw.WriteLine("1");
      sw.WriteLine("# Координаты Х");
      for (int i = 0; i < countCircle; i++)
      {
        sw.Write("{0} ", x[2 * i]);
      }
      sw.WriteLine();
      sw.WriteLine("# Координаты Y");
      for (int i = 0; i < countCircle; i++)
      {
        sw.Write("{0} ", x[2 * i + 1]);
      }
      sw.WriteLine();
      sw.WriteLine("# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ");
      sw.WriteLine("# Радиуc внешнего круга");
      sw.WriteLine("{0} ", x[2 * countCircle]);
      sw.Close();
    }

    static void SaveToFile(string fileName, double[] obj, int sizeC)
    {
      using (var sw = new StreamWriter(fileName, false))
      {
        sw.WriteLine("{0} ", obj[2 * sizeC]); // Ro
        sw.WriteLine("{0} ", sizeC);          //запись размерностей масива
        for (int i = 0; i < sizeC; i++)
          sw.WriteLine("{0} {1} {2}", obj[i], obj[2 * i], obj[2 * i + 1]);
      }
      //for (int i = 0; i < x.Length; ++i)
      //Console.WriteLine("x[{0}]={1}", i, x[i]);
    }
    static void OpenFromFile(string fileName, out double[] r, out int[] group, out Restriction rType, out double[] x, out double[] y)
    {
      bool tryXY = true;
      int cicleCount = 0;
      r = null;
      group = null;
      x = null;
      y = null;
      rType = Restriction.Power;
      // Create an instance of StreamReader to read from a file.
      // The using statement also closes the StreamReader.
      using (var sr = new StreamReader(fileName, false))
      {
        string line;
        try
        {
          // Блок: Тип задачи
          while ((line = sr.ReadLine()) != null)
          {
            if (line[0] == '#') continue;
            switch (Convert.ToInt32(line))
            {
              case 0:
                rType = Restriction.FixedRadius;
                break;
              case 1:
                rType = Restriction.Power;
                break;
              case 2:
                rType = Restriction.Lineal;
                break;
            }
            break;
          }
          // Блок: Количество кругов
          while ((line = sr.ReadLine()) != null)
          {
            if (line[0] == '#') continue;
            cicleCount = Convert.ToInt32(line);
            break;
          }
          r = new double[cicleCount];
          group = new int[cicleCount];
          // Блок: Радиусы
          while ((line = sr.ReadLine()) != null)
          {
            if (line[0] == '#') continue;
            line = line.Trim(' ');
            var elements = line.Split(new Char[] { ' ' });
            for (int i = 0; i < cicleCount; ++i)
              r[i] = Convert.ToDouble(elements[i]);
            break;
          }
          // Блок: Группы
          while ((line = sr.ReadLine()) != null)
          {
            if (line[0] == '#') continue;
            line = line.Trim(' ');
            var elements = line.Split(new Char[] { ' ' });
            for (int i = 0; i < cicleCount; ++i)
              group[i] = Convert.ToInt32(elements[i]);
            break;
          }
          // Блок: Будут ли задаваться координаты х и у
          while ((line = sr.ReadLine()) != null)
          {
            if (line[0] == '#') continue;
            line = line.Trim(' ');
            tryXY = (Convert.ToInt32(line) == 1) ? true : false;
            break;
          }
          if (tryXY)
          {
            // Блок: X
            while ((line = sr.ReadLine()) != null)
            {
              if (line[0] == '#') continue;
              line = line.Trim(' ');
              var elements = line.Split(new Char[] { ' ' });
              for (int i = 0; i < cicleCount; ++i)
                x[i] = Convert.ToDouble(elements[i]);
              break;
            }
            // Блок: Y
            while ((line = sr.ReadLine()) != null)
            {
              if (line[0] == '#') continue;
              line = line.Trim(' ');
              var elements = line.Split(new Char[] { ' ' });
              for (int i = 0; i < cicleCount; ++i)
                y[i] = Convert.ToDouble(elements[i]);
              break;
            }
          }
        }
        catch (Exception ex) { Console.WriteLine(ex.Message); }
        finally
        {
          sr.Dispose();
        }
      }
    }
    /// <summary>
    /// Отобразить входные данные
    /// </summary>
    /// <param name="data">Данные передаваемы в решатель: переменные</param>
    static void ShowData(double[] data)
    {
      int cicleCount = (data.Length - 1) / 3;
      Console.WriteLine("~~~~~~~~~~~ Координата Х ~~~~~~~~~~~");
      for (int i = 0; i < cicleCount; ++i)
        Console.WriteLine(" x[{0}]= {1}", i + 1, data[2 * i]);
      Console.WriteLine("~~~~~~~~~~~ Координата Y ~~~~~~~~~~~");
      for (int i = 0; i < cicleCount; ++i)
        Console.WriteLine(" y[{0}]= {1}", i + 1, data[2 * i + 1]);
      Console.WriteLine("~~~~~~~~~~~ Radius ~~~~~~~~~~~");
      for (int i = 0; i < cicleCount; ++i)
        Console.WriteLine(" r[{0}]= {1}", i + 1, data[2 * cicleCount + i]);
      Console.WriteLine("\n R = {0}", data[3 * cicleCount]);
    }
    /// <summary>
    /// Для задачи с фиксированным радиусами
    /// </summary>
    /// <param name="data">Данные передаваемы в решатель: переменные</param>
    /// <param name="radius">Константные радиусы</param>
    static void ShowData(double[] data, double[] radius)
    {
      int cicleCount = radius.Length;
      Console.WriteLine("~~~~~~~~~~~ Координата Х ~~~~~~~~~~~");
      for (int i = 0; i < cicleCount; ++i)
        Console.WriteLine(" x[{0}]= {1}", i + 1, data[2 * i]);
      Console.WriteLine("~~~~~~~~~~~ Координата Y ~~~~~~~~~~~");
      for (int i = 0; i < cicleCount; ++i)
        Console.WriteLine(" y[{0}]= {1}", i + 1, data[2 * i + 1]);
      Console.WriteLine("~~~~~~~~~~~ Radius ~~~~~~~~~~~");
      for (int i = 0; i < cicleCount; ++i)
        Console.WriteLine(" r[{0}]= {1}", i + 1, radius[i]);
      Console.WriteLine("\n R = {0}", data[2 * cicleCount]);
    }
  }
}
