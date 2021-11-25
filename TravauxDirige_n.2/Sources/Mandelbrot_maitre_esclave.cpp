#include <iostream>
#include <cstdlib>
#include <string>
#include <chrono>
#include <cmath>
#include <vector>
#include <fstream>
#include "mpi.h"



struct Complex
{
    Complex() : real(0.), imag(0.)
    {}
    Complex(double r, double i) : real(r), imag(i)
    {}
    Complex operator + (const Complex& z)
    {
        return Complex(real + z.real, imag + z.imag);
    }
    Complex operator * (const Complex& z)
    {
        return Complex(real * z.real - imag * z.imag, real * z.imag + imag * z.real);
    }
    double sqNorm() { return real * real + imag * imag; }
    double real, imag;
};

std::ostream& operator << (std::ostream& out, const Complex& c)
{
    out << "(" << c.real << "," << c.imag << ")" << std::endl;
    return out;
}

int iterMandelbrot(int maxIter, const Complex& c)
{
    Complex z{ 0.,0. };

    if (c.real * c.real + c.imag * c.imag < 0.0625)
        return maxIter;
    if ((c.real + 1) * (c.real + 1) + c.imag * c.imag < 0.0625)
        return maxIter;
   
    if ((c.real > -0.75) && (c.real < 0.5)) {
        Complex ct{ c.real - 0.25,c.imag };
        double ctnrm2 = sqrt(ct.sqNorm());
        if (ctnrm2 < 0.5 * (1 - ct.real / ctnrm2)) return maxIter;
    }
    int niter = 0;
    while ((z.sqNorm() < 4.) && (niter < maxIter))
    {
        z = z * z + c;
        ++niter;
    }
    return niter;
}


void
computeMandelbrotSetRow(int W, int H, int maxIter, int num_ligne, int* pixels)
{
    // Calcul le facteur d'échelle pour rester dans le disque de rayon 2
    // centré en (0,0)
    double scaleX = 3. / (W - 1);
    double scaleY = 2.25 / (H - 1.);
    //
    // On parcourt les pixels de l'espace image :
    for (int j = 0; j < W; ++j) {
        Complex c{ -2. + j * scaleX,-1.125 + num_ligne * scaleY };
        pixels[j] = iterMandelbrot(maxIter, c);
    }
}

std::vector<int>
computeMandelbrotSet(int W, int H, int maxIter, int nbp, int rank)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::vector<int> pixels(W * H);
    std::vector<int> pixels_local(W * int(H / nbp));
    start = std::chrono::system_clock::now();
;

    MPI_Status Stat, Stat2;

    if (rank == 0) {

        int count_task = 0;
        int np_tasks = H;
        for (int i = 1; i < nbp; i++) {
            MPI_Send(&count_task, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            count_task += 1;
        }

        
        int task_performed; 
        while (count_task < np_tasks) {
            MPI_Recv(&task_performed, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &Stat);
            MPI_Recv(pixels.data() + W * task_performed, W, MPI_INT, Stat.MPI_SOURCE, 0, MPI_COMM_WORLD, &Stat2); // rcv 2 msg de meme proc?
            MPI_Send(&count_task, 1, MPI_INT, Stat.MPI_SOURCE, 0, MPI_COMM_WORLD);
            count_task += 1;
        }

        count_task = -1;
        for (int i = 1; i < nbp; i++) {
            MPI_Send(&count_task, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        end = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "Temps calcul ensemble mandelbrot : " << elapsed_seconds.count()
            << std::endl;
    }
    else {
        MPI_Status Stat3;
        std::vector<int> pixels_line(W);
        int num_task = 0;
        while (num_task != -1) {
            MPI_Recv(&num_task, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Stat3);
            if(num_task >= 0){
                computeMandelbrotSetRow(W, H, maxIter, num_task, pixels_line.data());
                MPI_Send(&num_task, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                MPI_Send(pixels_line.data(), W, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

    return pixels;
}


void savePicture(const std::string& filename, int W, int H, const std::vector<int>& nbIters, int maxIter)
{
    double scaleCol = 1. / maxIter;
    std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
    ofs << "P6\n"
        << W << " " << H << "\n255\n";
    for (int i = 0; i < W * H; ++i) {
        double iter = scaleCol * nbIters[i];
        unsigned char r = (unsigned char)(256 - (unsigned(iter * 256.) & 0xFF));
        unsigned char b = (unsigned char)(256 - (unsigned(iter * 65536) & 0xFF));
        unsigned char g = (unsigned char)(256 - (unsigned(iter * 16777216) & 0xFF));
        ofs << r << g << b;
    }
    ofs.close();
}

int main(int argc, char* argv[])
{
    int nbp, rank;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nbp);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int W = 800;
    const int H = 600;

    const int maxIter = 8 * 65536;
    auto iters = computeMandelbrotSet(W, H, maxIter, nbp, rank);
    if(rank == 0){
        savePicture("mandelbrot.tga", W, H, iters, maxIter);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
