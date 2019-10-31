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
#include <functional>

using namespace std;


class injectedVect{
   public:
	vector<double> vect;
	vector<int> potential_injection_sites;
	
	injectedVect(vector<double> init_vect){
		vect = init_vect;
	}

	void add_injection_site(int index){
		potential_injection_sites.push_back(index);
	}

	void clear_injection_sites(){
		potential_injection_sites.clear();
	}

	int inject_failure(){
		if (potential_injection_sites.size() == 0)
			return -1;

		// randomly select a position and flip a random bit
		random_device random_device;
		mt19937 engine{random_device()};
        	uniform_int_distribution<int> dist(
				0, potential_injection_sites.size() - 1);
		int index = potential_injection_sites[dist(engine)];
		ieee754_double random_element = {vect[index]};
		random_element.ieee.mantissa1 |= 1u << 16; // set bit 16 of mantissa
		vect[index] = random_element.d;
		return index;
	}

	int get_total_corruptions(vector<double> correct_vector){
		unsigned int i, cnt = 0;
		for (i=0; i<correct_vector.size(); i++)
			if (correct_vector[i] != vect[i])
				cnt ++;
		return cnt;
	}
};

class CSRMatrix{
   public:
        int rows;
        int cols;
        int nnz;
        vector<double> values;
        vector<int> aj, ai;

	CSRMatrix(const char* file_path){
		string fn = file_path;
		if(fn.substr(fn.find_last_of(".") + 1) == "csr")
			read_matrix_from_csr(file_path);
		else
                	read_matrix_from_file(file_path);
	}

	void print(){
		cout << "Matrix stored in CSR format" << endl;
		cout << "Matrix non zero values" << endl;
		for (vector<double>::const_iterator i = values.begin(); i != values.end(); ++i)
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

        vector<double> multiply(vector<double> v){
		int vsize = v.size();
		assert(cols == vsize);

                vector<double> result(v.size());
		int r_index;
                for(int i = 0; i < rows; i++){
                        for(r_index = ai[i]; r_index < ai[i + 1]; r_index++){
				result[i] += values[r_index] * v[aj[r_index]];
			}
		}
		return result;
        }

	void read_matrix_from_csr(const char *file_path)
	{
                std::ifstream fin(file_path);
		
                fin >> cols >> rows >> nnz;
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
			values.push_back(data);
		}
	  	fin.close();
	}//read_csr

        void read_matrix_from_file(const char* file_path)
        {
                std::ifstream fin(file_path);
		
                fin >> rows >> cols >> nnz;
                ai.push_back(0);

                int last_row = -1;
                for (int i = 0; i < nnz; i++)
                {
                        int r, c;
                        double data;
                        fin >> r >> c >> data;
                        values.push_back(data);
                        aj.push_back(c);
                        if (last_row != r) 
				for (int j = 0; j < r-last_row; j++)
                                	ai.push_back(ai.back());
                        last_row = r;
                        ai.back() += 1;
		}

                fin.close();
        }//read_matrix
};

class spMV{
   public:
        int max_iterations;
	int step;
        vector<double> vect;
	CSRMatrix matrix;
	injectedVect faulty_vect;
	function<bool(double)> inject_criteria;

        spMV(vector<double> x,
	     const char* file_path,
	     const function<bool(double)> &fct,
	     int max_iter): matrix(file_path), faulty_vect(x), inject_criteria(fct){
                vect = x;
		step = 0;
                max_iterations = max_iter;
        }

	void print_vect(){
		for (vector<double>::const_iterator i = vect.begin(); i != vect.end(); ++i)
			cout << *i << ' ';
		cout << endl;
	}

	void print_matrix(){
		matrix.print();
	}

	int get_max_velocity(vector<double> v2){
		// identify vector elements of high changing velocity
		int max_velocity = 0;
		for(size_t i = 0; i < vect.size(); i++)
		{
			if(max_velocity < abs(v2[i] - vect[i]))
				max_velocity = abs(v2[i] - vect[i]);
		}
		return max_velocity;
	}

	int inject_failures(vector<double> v2){
		int max_velocity = get_max_velocity(v2);
		// if the elements did not changed
		if (max_velocity == 0)
			return -1;
		for(size_t i = 0; i < vect.size(); i++)
		{
			double velocity = abs(v2[i] - vect[i])/max_velocity;
			if (inject_criteria(velocity))
				faulty_vect.add_injection_site(i);
		}
		return faulty_vect.inject_failure();
	}
	

	int run(){
		vector<double> result;
		int i;
		for (i = 0; i < max_iterations; i++)
		{
			result = matrix.multiply(vect);
			faulty_vect = matrix.multiply(faulty_vect.vect);
			// if the solution converged
			if (get_max_velocity(result) == 0)
				break;
			step ++;
			if (step == 1)
				inject_failures(result);
			vect = result;
			cout << step << " failures: " << faulty_vect.get_total_corruptions(vect) << endl;
		}
		cout << "Final failure count: " << faulty_vect.get_total_corruptions(vect) << endl;
		return step;
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
        	vect.push_back(rand() % 10);
		cout << vect.back() << " ";
	}

	spMV sim(vect, argv[1], [](double val) { return (val < 0.7); }, 10);
	sim.print_matrix();

	int sim_steps = sim.run();
	if (sim_steps < 10)
		cout << "Converge in " << sim_steps << " steps" << endl;
	cout << "Multiplication result:"<<endl; 
	sim.print_vect();
}
