// Produit matrice-vecteur
# include <cassert>
# include <vector>
# include <iostream>
# include <mpi.h>

// ---------------------------------------------------------------------
class Matrix : public std::vector<double>
{
public:
    Matrix (int dim);
    Matrix( int nrows, int ncols );
    Matrix( const Matrix& A ) = delete;
    Matrix( Matrix&& A ) = default;
    ~Matrix() = default;

    Matrix& operator = ( const Matrix& A ) = delete;
    Matrix& operator = ( Matrix&& A ) = default;
    
    double& operator () ( int i, int j ) {
        return m_arr_coefs[i + j*m_nrows];
    }
    double  operator () ( int i, int j ) const {
        return m_arr_coefs[i + j*m_nrows];
    }
    
    std::vector<double> operator * ( const std::vector<double>& u ) const;
    
    std::ostream& print( std::ostream& out ) const
    {
        const Matrix& A = *this;
        out << "[\n";
        for ( int i = 0; i < m_nrows; ++i ) {
            out << " [ ";
            for ( int j = 0; j < m_ncols; ++j ) {
                out << A(i,j) << " ";
            }
            out << " ]\n";
        }
        out << "]";
        return out;
    }
private:
    int m_nrows, m_ncols;
    std::vector<double> m_arr_coefs;
};
// ---------------------------------------------------------------------
inline std::ostream& 
operator << ( std::ostream& out, const Matrix& A )
{
    return A.print(out);
}
// ---------------------------------------------------------------------
inline std::ostream&
operator << ( std::ostream& out, const std::vector<double>& u )
{
    out << "[ ";
    for ( const auto& x : u )
        out << x << " ";
    out << " ]";
    return out;
}
// ---------------------------------------------------------------------
std::vector<double> 
Matrix::operator * ( const std::vector<double>& u ) const
{
    MPI_Comm globComm;
    MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
    int nbp;
    MPI_Comm_size(globComm, &nbp);
    int rank;
    MPI_Comm_rank(globComm, &rank);
    int pas = int (m_nrows/nbp);

    const Matrix& A = *this;
    assert( u.size() == unsigned(m_ncols) );
    std::vector<double> v(pas, 0.);
    for ( int i = rank*pas; i < (rank+1)*pas; ++i ) 
    {
        for ( int j = 0; j < m_ncols; ++j ) 
        {
            v[i - rank*pas] += A(i,j)*u[j];
        }            
    }


    MPI_Barrier(MPI_COMM_WORLD);
    double recbuf [m_nrows];
    MPI_Gather (v.data (), pas, MPI_DOUBLE, recbuf, pas, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    std::cout << std::endl << "here" << v << std::endl << std::endl;
    return (vector<double> (recbuf, recbuf+m_nrows));
    
    
    
}

// =====================================================================
Matrix::Matrix (int dim) : m_nrows(dim), m_ncols(dim),
                           m_arr_coefs(dim*dim)
{
    for ( int i = 0; i < dim; ++ i ) {
        for ( int j = 0; j < dim; ++j ) {
            (*this)(i,j) = (i+j)%dim;
        }
    }
}
// ---------------------------------------------------------------------
Matrix::Matrix( int nrows, int ncols ) : m_nrows(nrows), m_ncols(ncols),
                                         m_arr_coefs(nrows*ncols)
{
    int dim = (nrows > ncols ? nrows : ncols );
    for ( int i = 0; i < nrows; ++ i ) {
        for ( int j = 0; j < ncols; ++j ) {
            (*this)(i,j) = (i+j)%dim;
        }
    }    
}
// =====================================================================
int main( int nargs, char* argv[] )
{
    MPI_Init( &nargs, &argv );
    const int N = 120;
    Matrix A(N);
    int rang;
    MPI_Comm_rank(MPI_COMM_WORLD, &rang);
    if (rang == 0)
    {
        std::cout  << "A : " << A << std::endl;        
    }
    std::vector<double> u( N );
    for ( int i = 0; i < N; ++i ) u[i] = i+1;
    if (rang == 0)
        std::cout << " u : " << u << std::endl;
    std::vector<double> v = A*u;
    if (rang == 0)
        std::cout << "process :" << rang << "A.u = " << v << std::endl;

    MPI_Finalize ();
    return EXIT_SUCCESS;
}
