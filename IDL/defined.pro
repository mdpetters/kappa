; defined: this function is just uses the IDL function size to determine
; if a variable is defined or not
function defined,variable

  test = size(variable)
  if test(n_elements(test )-2) NE 0 then return,1
  return,0
end
