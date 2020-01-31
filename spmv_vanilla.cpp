#include <iostream>
#include <fstream>
#include <vector> 
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <sys/time.h>
#include <cmath>
#include <random>
#include <bitset> 
#include <functional>
#include <algorithm>

using namespace std;

static_assert(sizeof(double) == sizeof(unsigned long long), "");

union Converter { uint64_t i; double d; };

class injectedVect{
   public:
	vector<double> vect;
	vector<int> potential_injection_sites;
	size_t verbose = 0;
	
	injectedVect(vector<double> init_vect){
		vect = init_vect;
	}

	void add_injection_site(int index){
		potential_injection_sites.push_back(index);
	}

	void clear_injection_sites(){
		potential_injection_sites.clear();
	}

	void print_vect(){
		for (vector<double>::const_iterator i = vect.begin(); i != vect.end(); ++i)
			cout << *i << ' ';
		cout << endl;
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
		
		Converter c;
		c.d = vect[index];
	    bitset<sizeof(double) * 8> b(c.i);
		if (verbose){
			cout << "Double value  : " << vect[index] << endl;
			cout << "Int value : " << b.to_ulong() << endl;
			cout << "BitSet : " << b.to_string() << endl;
		}
	   	
		b.flip(57); // flip the middle exponent bit
		//b.flip(52); // flip the most insignificant exponent bit
		c.i = b.to_ullong();
		
		if (verbose){ 
			cout << "BitSet : " << b.to_string() << endl;
	    	cout << b.to_ulong() << endl;
	    	cout << c.d << endl;
		}	
		
		vect[index] = c.d;
		return index;
	}

	int get_total_corruptions(vector<double> correct_vector){
		unsigned int i, cnt = 0;
		for (i=0; i<correct_vector.size(); i++)
			if (abs(correct_vector[i] - vect[i]) > 0.001){
				cnt ++;
			}
		return cnt;
	}

	double get_max_velocity(vector<double> v2){
		double max_velocity = 0;
		for(size_t i = 0; i < vect.size(); i++)
		{
			if(max_velocity < abs(v2[i] - vect[i]))
				max_velocity = abs(v2[i] - vect[i]);
		}
		return max_velocity;
	}
};

class resilientGradient{
	private:
		int iteration = 0;
		vector<double> previous1;
		vector<double> previous2;
	public:
		resilientGradient(vector<double> v): previous1(v) {}
		
		resilientGradient(vector<double> v1,
						  vector<double> v2): previous1(v1), previous2(v2) {}
	
		vector<double> first_gradient(vector<double> current, int iter){
			assert(current.size() == previous1.size());
			assert(iter > iteration);
			vector<double> res;
			for(size_t i = 0; i < current.size(); i++)
				res.push_back(abs(current[i] - previous1[i])/(iter-iteration));
			return res;
		}
	
		vector<double> second_gradient(vector<double> current, int iter){
			assert(current.size() == previous1.size());
			assert(iter > iteration);
			vector<double> res;
			for(size_t i = 0; i < current.size(); i++)
				res.push_back(abs(current[i] - previous1[i])/(iter-iteration));
			return res;
		}
};

class Vector{
	public:
		vector<double> vect;

		Vector(size_t size){
			randomize(size);
		}

		Vector(vector<double> v){
			vect = v;
		}
		
		void randomize(size_t size){
			vect.clear();
			for (size_t i =0; i < size; i++)
				vect.push_back(rand() % 100);
			normalize();
		}

		void update(vector<double> v){
				vect = v;
				normalize();
		}

		void normalize(){
			double sum = 0;
			size_t i;
			for (i=0; i<vect.size(); i++)
				sum += (vect[i] * vect[i]);
			sum = sqrt(sum);
			for (i=0; i<vect.size(); i++)
				vect[i] /= sum;
		}

		double norm_diffrence(vector<double> v){
			double sum = 0;
			size_t i;
			for (i=0; i<vect.size(); i++)
				sum += ((v[i] - vect[i]) * (v[i] - vect[i]));
			sum = sqrt(sum);
			return sum;
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
		double sum = 0;
        for(int i = 0; i < rows; i++){
            for(r_index = ai[i]; r_index < ai[i + 1]; r_index++){
				result[i] += values[r_index] * v[aj[r_index]];
			}
			sum += (result[i] * result[i]);
		}
		sum = sqrt(sum);
		for (int i =0; i < rows; i++)
			result[i] /= sum;

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
			ai.push_back(data - 1);
	    }
		for(i=0; i < nnz; ++i)
		{
            int data;
            fin >> data;
			aj.push_back(data - 1);
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
                aj.push_back(c - 1);
				r = r - 1;
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
	private:
		vector<int> failures_per_step;
		vector<int> const_elements_per_step;
		vector<int> slow_elements_per_step;
	public:
		int max_iterations;
		int step;
		vector<double> vect;
		CSRMatrix matrix;
		injectedVect faulty_vect;
		function<bool(double)> inject_criteria;
		string file_path;

        spMV(vector<double> x, const char* file_path,
	     const function<bool(double)> &fct, int max_iter,
             string fp): matrix(file_path), faulty_vect(x), inject_criteria(fct), file_path(fp){
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

	void write_failure_to_file(int velocity){
		ofstream fout(file_path + ".failures", ios::app);
		fout << velocity << ' ';
		
		for (vector<int>::const_iterator i = failures_per_step.begin(); i != failures_per_step.end(); ++i)
			fout << *i << ' ';
		fout << endl;
	}

	void write_slow_elements_to_file(int velocity){
		ofstream fout_c(file_path + ".constant", ios::app);
		ofstream fout_s(file_path + ".slow", ios::app);
		fout_c << velocity << ' ';
		fout_s << velocity << ' ';
		for (vector<int>::const_iterator i = const_elements_per_step.begin(); i !=const_elements_per_step.end(); ++i)
			fout_c << *i << ' ';
		for (vector<int>::const_iterator i = slow_elements_per_step.begin(); i !=slow_elements_per_step.end(); ++i)
			fout_s << *i << ' ';
		fout_c << endl;
		fout_s << endl;
	}

	double get_max_velocity(vector<double> v2){
		// identify vector elements of high changing velocity
		double max_velocity = 0;
		for(size_t i = 0; i < vect.size(); i++)
		{
			if(max_velocity < abs(v2[i] - vect[i]))
				max_velocity = abs(v2[i] - vect[i]);
		}
		return max_velocity;
	}

	int inject_failures(vector<double> v2){
		double max_velocity = get_max_velocity(v2);
		// if the elements did not changed
		if (max_velocity < 0.001)
			return -1;
		for(size_t i = 0; i < vect.size(); i++)
		{
			double velocity = abs(v2[i] - vect[i])/max_velocity; 
			if (inject_criteria(velocity)){
				faulty_vect.add_injection_site(i);
			}
		}
		return faulty_vect.inject_failure();
	}
	

	void store_const_elem(vector<double> v){
		int con=0, slow=0;
		for (size_t i=0; i< v.size(); i++){
			if (abs(v[i]-vect[i])<0.001){
				con ++;
				continue;
			}
			if (abs(v[i]-vect[i])<0.01){
				slow++;
			}
		}
		const_elements_per_step.push_back(con);
		slow_elements_per_step.push_back(slow);
	}

	int run(){
		vector<double> result;
		int i;
		for (i = 0; i < max_iterations; i++)
		{
			// multiply the faulty vector
			result = matrix.multiply(faulty_vect.vect);
			// if the solution converged
			if (faulty_vect.get_max_velocity(result) < 0.001)
				break;
			faulty_vect.vect = result;
			// multiply the correct matrix
			result = matrix.multiply(vect);
			step ++;
			if (step == 2)
				inject_failures(result);
			store_const_elem(result);
			vect = result;

			failures_per_step.push_back(faulty_vect.get_total_corruptions(vect));
		}

		return step;
	}
};

function<bool(double)> get_injection_boundries(int param){
	switch(param){
		case 1:
			return [](double val) { return (val >= 0) && (val < 0.1); };
		case 2:
			return [](double val) { return (val >= 0.1) && (val < 0.5); };
		case 3:
			return [](double val) { return (val >= 0.5) && (val < 0.9); };
		case 4:
			return [](double val) { return (val >= 0.9) && (val <= 1); };
	}
	return [](double val) { return (val < 0); };
}

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		cout << "Usage: " << argv[0] << " matrix_input_file vector_size output_file" << endl;
		return 1;
	}
	bool verbose = 0;
	int v, size = strtol(argv[2], NULL, 10);
	if (verbose)
		cout << "Problem size: " << size << endl;

	srand((unsigned) time(NULL));
	Vector vect(size);
	for (int loop = 0; loop < 1; loop ++){
		// create random vector (values between 0 and 100)
		vect.randomize(size);

		// simulate injection for different velocities
		for (v=2; v<3; v++){
			function<bool(double)> injection_fct = get_injection_boundries(v);
			
			// create simulation environment
			spMV sim(vect.vect, argv[1], injection_fct, 1000, argv[3]);
			if (verbose)
				sim.print_matrix();

			int sim_steps = sim.run();

			// write the number of failed sites per steps
			sim.write_failure_to_file(v);
			sim.write_slow_elements_to_file(v);
			//cout << "Convergence in " << sim_steps << " steps" << endl;
			cout << v << ' ' << sim_steps << endl;
		
			if (verbose){
				cout << endl << "Multiplication result:"<<endl; 
				sim.print_vect();
			}
		}
	}
	cout << endl;
}
