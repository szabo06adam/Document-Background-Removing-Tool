#include <cstdlib>
#include <stdlib.h>
#include <iostream>

#include "Document.h"

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>


/* Usage example */

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
	/* Options: */
	bool gsc_output = true;
	std::string gsc_output_name = "output.png";
	bool rgb_output = false;
	std::string rgb_output_name = "";
	bool scale_rgb = true;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("input-file,I", po::value<std::string>(), "The input file.")
		("output-file,O", po::value<std::string>(), "The greyscale output file.")
		("no-greyscale,N", "No greyscale image output.")
		("colour,C", po::value<std::string>(), "Colour image output.")
		("dont-scale-colour,D", "Do not scale lighting on colour image.");

	po::positional_options_description pd;
	pd.add("input-file", 1);
	pd.add("output-file", 1);

	/* Parsing: */
	po::variables_map vm;
	po::store(po::command_line_parser(argc,argv).options(desc).positional(pd).run(), vm);
	po::notify(vm);    

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 1;
	}

	if(!vm.count("input-file"))
	{
		std::cout << "No input image file." << std::endl;
		exit(EXIT_SUCCESS);
	}
	if(vm.count("output-file"))
		gsc_output_name = vm.at("output-file").as<std::string>();

	if (vm.count("colour"))
	{
		rgb_output_name = vm["colour"].as<std::string>();
		rgb_output = true;
		if(vm.count("no-greyscale"))
			gsc_output = false;
		if(vm.count("dont-scale-colour"))
			scale_rgb = false;
	}
	else
	{
		if(vm.count("no-greyscale"))
			std::cout << "--no-greyscale option can only be used with -colour." << "\n";
		if(vm.count("dont-scale-colour"))
			std::cout << "--dont-scale-colour option can only be used with -colour." << "\n";
	}

	Document::init();
	Document doc = Document(vm["input-file"].as<std::string>());
	doc.RGBtoGSc();
	doc.InitZones();
	doc.ScaleLighting(scale_rgb);
	doc.DeleteBackground();
	if(rgb_output) doc.WriteRGB(rgb_output_name);
	if(gsc_output) doc.WriteGSc(gsc_output_name);

	return(EXIT_SUCCESS);
}