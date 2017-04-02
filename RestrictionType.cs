using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace hs071_cs
{
  // Классы задач 
  public enum Restriction
  {
    Power /*степенные ограничения*/,
    Lineal /*линейные ограничения*/,
    AlphaPower,
    FixedRadius
  }
  public enum CircleTaskTypes
  {
    FixedRadius,
    DynamicRadius
  }
}
