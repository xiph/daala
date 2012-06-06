#if !defined(_cpack_linalg_svd_H)
# define _cpack_linalg_svd_H (1)

/*Computes the singular value decomposition of a matrix into matrices
   U\Sigma V^T, where U and V are orthonormal and \Sigma is diagonal.
  _w:     An _nrows by (_nrows+_ncols) matrix of storage for input and output.
          On input, the first _nrows rows contain the input matrix.
          On output, the first _nrows rows contain the columns of U, each
           scaled by the corresponding singular value.
          The caller must divide by the appropriate value to obtain a unit
           vector, if one is desired.
          The next _ncols rows contain the columns of V (_not_ V^T).
  _s:     On output contains the squares of the _ncols singular values.
  _nrows: The number of rows of the matrix.
  _ncols: The number of columns of the matrix.
  Return: The estimated column rank of the matrix.*/
int svd(double **_w,double *_s,int _nrows,int _ncols);
/*Computes the Moore-Penrose pseudoinverse of a matrix using an SVD.
  _w:     An _nrows by (_nrows+_ncols) matrix of storage for input and output.
          On input, the first _nrows rows contain the input matrix.
          On output, the first _nrows rows contain the transpose of the
           pseudoinverse.
          The next _ncols rows are temporary storage.
          On output the contents are undefined.
  _s:     _ncols temporary values.
          On output the contents are undefined.
  _nrows: The number of rows of the matrix.
  _ncols: The number of columns of the matrix.
  Return: The estimated column rank of the matrix.*/
int svd_pseudoinverse(double **_w,double *_s,int _nrows,int _ncols);

#endif                                                   /*_cpack_linalg_svd_H*/
