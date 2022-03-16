def getValue(value, vtype):
  if vtype == "double":
      return "%.15e"%value
  elif vtype == "float":
      return "%2.10e"%value
  else:
      return value
