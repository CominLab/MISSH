//============================================================================
// Name        : MISSH.cpp
//============================================================================

#include "Test/Test.h"
#include "cxxopts.hpp"


int main(int argc, char* argv[]) {
	string dir_output = "../output/";
	bool sequence = false;
	FileParameter param;
	omp_set_num_threads(4);	
	int test_kind = 2;


	
	cxxopts::Options options(argv[0], "Efficient Hashing of Multiple Spaced Seeds with Application");

	options.add_options()
		("s,si", "Input filename single-end", cxxopts::value<std::string>())
		("p,pi", "Input filenames paired-end", cxxopts::value<std::vector<std::string>>())
		("q,ss", "Spaced seeds path", cxxopts::value<std::string>())
		("d,dirO", "Output directory", cxxopts::value<std::string>())
		("n,num", "Number of previous hashes", cxxopts::value<int>())
		("t,test", "Test kind (single or multi)", cxxopts::value<std::string>())
		("m,threads", "Number of threads", cxxopts::value<int>())
		("h,help", "Print help");

	
	try
	{
		auto result = options.parse(argc, argv);

		if (result.count("help")) {
			std::cout << options.help() << std::endl;
			return 0;
		}

		if (result.count("si")) {
			std::string input = result["si"].as<std::string>();
			if (!param.init(input, "")) {
				std::cerr << "Please enter an input filename single-end: -si <AbsPathFile>" << std::endl;
				std::cout << options.help() << std::endl;
				return 0;
			}
			sequence = true;
		}

		if (result.count("pi")) {
			std::vector<std::string> inputs = result["pi"].as<std::vector<std::string>>();
			if (inputs.size() != 2 || !param.init(inputs[0], inputs[1])) {
				std::cerr << "Please enter input filenames paired-end: -pi <AbsPathFile1>,<AbsPathFile2>" << std::endl;
				std::cout << options.help() << std::endl;
				return 0;
			}
			sequence = true;
		}

		// Con questo comando è possibile fare in modo che invece di recuperare
		// tutte le posizioni da hash precedenti l'algoritmo ISSH recuperi da un
		// numero definito di hash precedenti il numero massimo di posizioni
		// recuperabili con quegli hash.
		// Attenzione! Deve essere impostato prima di -q, altrimenti non ha alcun
		// effetto.
		if (result.count("num")) {
			int numPrev = result["num"].as<int>();
			if (numPrev <= 0) {
				std::cerr << "Please enter valid number of previous hashes from which to retrieve positions." << std::endl;
				std::cout << options.help() << std::endl;
				return 0;
			} else {
				param.setNumPrev(static_cast<size_t>(numPrev));
			}
		}

		if (result.count("q")) {
			std::string pathQmers = result["q"].as<std::string>();
			std::vector<std::string> lines;
			if (getLines(pathQmers, lines)) {
				std::vector<bool> correctQmer(lines.size(), false);
				std::regex rgx("^1(0|1)*1$");
				for (size_t j = 0; j < lines.size(); j++) {
					correctQmer[j] = std::regex_match(lines[j], rgx);
					if (!correctQmer[j]) {
						std::cerr << "Error on " << j + 1 << "° spaced seed. Enter q-mer with 1 at begin and end of the string on input files. "
								<< "Ex. 1**1*11*1. 1 is the symbol considered, any others are not valid symbols." << std::endl;
						return 0;
					} else {
						param.addSpacedQmer(lines[j], lines[j]);
					}
				}
			} else {
				std::cerr << "Please enter a spaced seeds path as -q <AbsPathFile>. Every file's line must contain a spaced seeds." << std::endl;
				std::cout << options.help() << std::endl;
				return 0;
			}
		}

		if (result.count("dirO")) {
			dir_output = result["dirO"].as<std::string>();
			if (dir_output.empty()) {
				std::cerr << "Enter valid path for output directory." << std::endl;
				std::cout << options.help() << std::endl;
				return 0;
			}
		}

		if (result.count("test")) {
			std::string testValue = result["test"].as<std::string>();
			if (testValue == "single") {
				test_kind = 0;
			} else if (testValue == "multi") {
				test_kind = 1;
			} else {
				std::cerr << "Please enter valid value (single or multi). Do not write -test and its option to perform both tests." << std::endl;
				std::cout << options.help() << std::endl;
				return 0;
			}
		}

		if (result.count("threads")) {
			int threads_num = result["threads"].as<int>();
			if (threads_num <= 0) {
				std::cerr << "Please enter valid number of threads." << std::endl;
				std::cout << options.help() << std::endl;
				return 0;
			} else {
				omp_set_num_threads(threads_num);
			}
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error parsing options: " << e.what() << std::endl << std::endl;
		std::cout << options.help() << std::endl;
		return 1;
	}



	if(!sequence)
	{
		cerr<<endl<<"Please enter a DNA sequence with -si <AbsPathFile> or -pi <AbsPathFile1> <AbsPathFile2>\n"<<flush;
		return 0;
	}
	
	if(param.getVSpaced().empty())//Applied default
	{
		param.addSpacedQmer("CLARK-S", "1111011101110010111001011011111");
		param.addSpacedQmer("CLARK-S", "1111101011100101101110011011111");
		param.addSpacedQmer("CLARK-S", "1111101001110101101100111011111");
		param.addSpacedQmer("rasbhari_minimizing_overlap_complexity", "1111010111010011001110111110111");
		param.addSpacedQmer("rasbhari_minimizing_overlap_complexity", "1110111011101111010010110011111");
		param.addSpacedQmer("rasbhari_minimizing_overlap_complexity", "1111101001011100111110101101111");
		param.addSpacedQmer("rasbhari_maximizing_sensitivity", "1111011110011010111110101011011");
		param.addSpacedQmer("rasbhari_maximizing_sensitivity", "1110101011101100110100111111111");
		param.addSpacedQmer("rasbhari_maximizing_sensitivity", "1111110101101011100111011001111");
		cout << endl << "Applied default spaced seed" << flush;
	}

	//Creo cartella output se non presente
	createDirAndSubDir(dir_output);

	cout << endl << "Applied the following spaced seed..." << flush;
	for(size_t i = 0; i < param.getVSpaced().size(); i++)
		cout << endl << "Type:" << param.getVSpaced()[i].first << ", Spaced seed: " << param.getVSpaced()[i].second.toString() << flush;


	// Test for a single spaced seed at a time
	if (test_kind == 0 || test_kind == 2)
	{
		cout << endl << "Performing test for a single spaced seed at a time"  << flush;

		bool single_test_equals = false;
		string dir_output_1 = dir_output + "single/";
		Test test_single;
		if(test_single.load_sequences(param))
		{
			for(size_t j = 0; j < param.getVSpaced().size(); ++j)
			{
				string dir_output_tmp = dir_output_1 + param.getVSpaced()[j].first + param.getVSpaced()[j].second.toString() + "/";

				test_single.single_run(param.getVSpaced()[j].second, single_test_equals);
				test_single.single_save(param, dir_output_tmp);
			}
		}
	}

	// Test for a group of spaced seeds
	if (test_kind == 1 || test_kind == 2)
	{
		cout << endl << "Performing test for multiple spaced seeds at a time"  << flush;

		bool multi_test_equals = false;
		string dir_output_2 = dir_output + "multi/";
		// multi_spaced vector which contains all the spaced seeds given in input. 
		vector<SpacedQmer> multi_spaced;
		for(size_t i = 0; i < param.getVSpaced().size(); ++i)
			multi_spaced.push_back(param.getVSpaced()[i].second);

		Test test_multi;
		if(test_multi.load_sequences(param))
		{
			test_multi.multi_run(multi_spaced, multi_test_equals);
			test_multi.multi_save(param, multi_spaced, dir_output_2);
		}
	}

	cout<<"End\n";


	return 0;
}

