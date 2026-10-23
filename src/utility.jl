# 数值比较容差（Ref：可在运行时通过 set_quapi_tol! 调整）
const QuAPI_tol = Ref(1.0e-6)

"返回当前的数值比较容差。"
get_quapi_tol() = QuAPI_tol[]

"设置数值比较容差（返回新值）。"
set_quapi_tol!(x::Real) = (QuAPI_tol[] = float(x); QuAPI_tol[])
