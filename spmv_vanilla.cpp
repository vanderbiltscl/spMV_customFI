#include <iostream>
#include <fstream>
#include <vector> 
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <sys/time.h>
#include <cmath>
#include <random>
#include <ieee754.h>

using namespace std;

double timer() {
        struct timeval tp;
        gettimeofday(&tp, NULL);
        return ((double) (tp.tv_sec) + tp.tv_usec * 1e-6);
}

class injectedVect{
   public:
	vector<double> corrupted_vect;
	vector<int> potential_injection_sites;
	
	injectedVect(){
	}

	injectedVect(vector<double> init_vect){
		corrupted_vect = init_vect;
	}

	void set_init_vect(vector<double> init_vect){
		corrupted_vect = init_vect;
	}

	void add_injection_site(int index){
		potential_injection_sites.push_back(index);
	}

	void clear_injection_sites(){
		potential_injection_sites.clear();
	}

	int inject_failure(){
		// randomly select a position and flip a random bit
		random_device random_device;
		mt19937 engine{random_device()};
        	uniform_int_distribution<int> dist(
				0, potential_injection_sites.size() - 1);
	  	int index = potential_injection_sites[dist(engine)];
		ieee754_double random_element = {corrupted_vect[index]};
		random_element.ieee.mantissa1 |= 1u << 16; // set bit 16 of mantissa
		corrupted_vect[index] = random_element.d;
		return index;
	}

	int get_total_corruptions(vector<double> correct_vector){
		unsigned int i, cnt = 0;
		for (i=0; i<correct_vector.size(); i++)
			if (correct_vector[i] != corrupted_vect[i])
				cnt ++;
		return cnt;
	}
};

class spMV{
   public:
        int max_iterations;
        int rows;
        int cols;
        int nnz;
	int step;
        vector<double> matrix;
        vector<int> aj, ai;
        vector<double> vect;

        spMV(vector<double> x, const char* file_path){
                cols = x.size();
                vect = x;
		step = 0;
                max_iterations = 100;
		string fn = file_path;
		if(fn.substr(fn.find_last_of(".") + 1) == "csr")
			read_matrix_from_csr(file_path);
		else
                	read_matrix_from_file(file_path);
        }

	void print_vector(){
		cout << "Vector at step" << step << endl;
		for (vector<double>::const_iterator i = vect.begin(); i != vect.end(); ++i)
			cout << *i << ' ';
		cout << endl;
	}

	void print_matrix(){
		cout << "Matrix stored in CSR format" << endl;
		cout << "Matrix non zero values" << endl;
		for (vector<double>::const_iterator i = matrix.begin(); i != matrix.end(); ++i)
			cout << *i << ' ';
		cout << endl;
		cout << "Non zero values columns" << endl;
		for (vector<int>::const_iterator i = aj.begin(); i != aj.end(); ++i)
			cout << *i << ' ';
		cout << endl;
		cout << "Row indeces" << endl;
		for (vector<int>::const_iterator i = ai.begin(); i != ai.end(); ++i)
			cout << *i << ' ';
		cout << endl;
	}

	void add_over_vector(vector<double> v2){
		assert (v2.size() == vect.size());
		for(int i = 0; i < rows; i++)
			vect[i] += v2[i];
	}

        vector<double> multiply(){
                vector<double> result(cols);
		int r_index;
                for(int i = 0; i < rows; i++){
                        for(r_index = ai[i]; r_index < ai[i + 1]; r_index++){
				result[i] += matrix[r_index] * vect[aj[r_index]];
			}
		}
		return result;
        }

	void read_matrix_from_csr(const char *file_path)
	{
                std::ifstream fin(file_path);
		
		int n_cols;
                fin >> n_cols >> rows >> nnz;
                assert(cols == n_cols); 
                ai.push_back(0);

		int i;
	  	for(i=0; i < cols+1; ++i)
	    	{
                        int data;
                        fin >> data;
			ai.push_back(data);
	    	}
		for(i=0; i < nnz; ++i)
		{
                        int data;
                        fin >> data;
			aj.push_back(data);
		}
		for(i=0; i < nnz; ++i)
		{
                        double data;
                        fin >> data;
			matrix.push_back(data);
		}
	  	fin.close();
	}//read_csr

        void read_matrix_from_file(const char* file_path)
        {
                std::ifstream fin(file_path);
		
		int n_cols;
                fin >> rows >> n_cols >> nnz;
                assert(cols == n_cols); 
                ai.push_back(0);

                int last_row = -1;
                for (int i = 0; i < nnz; i++)
                {
                        int r, c;
                        double data;
                        fin >> r >> c >> data;
                        matrix.push_back(data);
                        aj.push_back(c);
                        if (last_row != r) 
				for (int j = 0; j < r-last_row; j++)
                                	ai.push_back(ai.back());
                        last_row = r;
                        ai.back() += 1;
                }

                fin.close();
        }//read_matrix

	vector<double> run(){
		vector<double> result;
		int i;
		for (i = 0; i < max_iterations; i++)
		{
			result = multiply();
			add_over_vector(result);
		}
		return vect;
	}
};

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		cout << "Usage: " << argv[0] << " matrix_file vector_size" << endl;
		return 1;
	}
	vector<double> vect;
	int size = strtol(argv[2], NULL, 10);
	cout << "Size: " << size << endl;
	srand((unsigned) time(NULL));
	for (int i =0; i < size; i++){
        	vect.push_back(1);//rand() % 10);
		cout << vect.back() << " ";
	}

	spMV sim(vect, argv[1]);
	sim.print_matrix();

	vector<double> result;
	result = sim.run();
	cout << "Multiplication result:"; 
	for (vector<double>::const_iterator i = result.begin(); i != result.end(); ++i)
		cout << *i << ' ';
	cout << endl;
}
