import config as cfg

class ScriptReader(object):
  @staticmethod
  def getColumns(filename):
    f = open(filename)
    result = []
    for value in f.readlines():
      result.append(value.strip())
    f.close()
    return result  
        
  @staticmethod
  def getHeader(filename):
    precision = cfg.floatPrecision
    f = open(filename)
    result = []
    for value in f.readlines():
        if value[0] == '#':
          if '#precision:' in value:
            tmp =  value.split(':')[1].strip()
            if tmp == cfg.doublePrecision :
              precision = tmp
          else:
            result.append(value)
    f.close()
    if len(result) == 0:
      return precision, None
    return precision, result  
