Para el proceso de compilación y ejecución del código usando Arrayfire +
CUDA + Armadillo es importante crear un enlace simbólico de la librería
libnvvm.so.3 (/usr/local/cuda/nvvm/lib64/libnvvm.so.3) a /usr/local/lib
esto con el fin de que Arrayfire pueda utilizarla.
