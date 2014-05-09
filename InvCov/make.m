function make()

mex QUIC.cpp QUIC_utils.cpp libmwlapack.lib libmwblas.lib

mex ML_QUIC.cpp QUIC_utils.cpp libmwlapack.lib libmwblas.lib

end

