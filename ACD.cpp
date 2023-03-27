#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <thread>
#include <string>
#include <vector>
#include <cmath>

int num = 1;               //Test number
double wA, wB, wC, wD;     //Frequencies
double k1, k3, k1d, k2d;   //Kappas
double p1 = 0.287;         //period(mm) double
double p2 = k1d;           //K1D double (1334.3984)
double p3 = k2d;           //K2D double (386.8075)
double p4 = 38.0;          //Krad double
double p5 = 0.0;           //Ki double
double p6 = 0.0;           //reflectivity double
int p7 = 50;               //Time-step

//File storing path:
// *** Mac ***
// std::string p8 = "/Users/joshua/Desktop/ACD/PMC.csv";
// std::string indPath = "/Users/joshua/Desktop/ACD/designs/xxx.ind";
// std::string settingsPath = "/Users/joshua/Desktop/ACD/Settings.txt";
// std::string textBlockPath = "/Users/joshua/Desktop/ACD/TextBlock.txt";
// std::string pweDataPath = "/Users/joshua/Desktop/ACD/bin/bstmp_eigs_TM.dat";
// std::string exePath = "/Users/joshua/Desktop/PMC_MacOS.exe";
// std::string recordPath = "/Users/joshua/Desktop/ACD/data.csv";

// *** Win ***
//std::stringstream n << setw(4) << setfill('0') << num;
std::string p8 = "C:\\ACD\\Tests\\PMC.csv";
std::string indPath = "C:\\ACD\\designs\\Test01.ind"; //+ n.str()
std::string settingsPath = "C:\\ACD\\Settings.txt";
std::string textBlockPath = "C:\\ACD\\TextBlock.txt";
std::string pweDataPath = "C:\\ACD\\bin\\bstmp_eigs_TM.dat";
std::string exePath = "C:\\ACD\\bin\\PMC.exe";
std::string recordPath = "C:\\ACD\\data.csv";
std::string bandSolvePath = "C:\\Synopsys\\PhotonicSolutions\\2022.06\\RSoft\\bin\\bandsolve.exe";


/* ① */
void createIndFile(){   
    // Create i/o file stream objects and open the file
    std::ofstream outfile(indPath);
    std::ifstream settings(settingsPath);
    std::ifstream textBlock(textBlockPath);

    /* Read and write Settings.txt line by line into outfile */
    if (settings.is_open()) {
        if (outfile.is_open()) {
            std::string line;
            while (std::getline(settings, line)) {
                outfile << line << std::endl;
            }
        }
        else {
            std::cerr << "Error: Could not open output file." << std::endl;
        }
    }
    else {
        std::cerr << "Unable to open settings.txt." << std::endl;
    }
    settings.close();   

    /* Generate polygons */
    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);   

    // Generate random 4 coordinates
    double point1x = (dist(gen)-0.5)*0.4;
    double point1y = (dist(gen)-0.5)*0.4;

    double point2x = (dist(gen)-0.5)*0.4;
    double point2y = (dist(gen)-0.5)*0.4;

    double point3x = (dist(gen)-0.5)*0.4;
    double point3y = (dist(gen)-0.5)*0.4;

    double point4x = (dist(gen)-0.5)*0.4;
    double point4y = (dist(gen)-0.5)*0.4;

   std::stringstream p1x;
   std::stringstream p1y;

   std::stringstream p2x;
   std::stringstream p2y;

   std::stringstream p3x;
   std::stringstream p3y;

   std::stringstream p4x;
   std::stringstream p4y;

   //Set precision and then convert to stringstream
   p1x << std::fixed << std::setprecision(3) << point1x;
   p1y << std::fixed << std::setprecision(3) << point1y;

   p2x << std::fixed << std::setprecision(3) << point2x;
   p2y << std::fixed << std::setprecision(3) << point2y;

   p3x << std::fixed << std::setprecision(3) << point3x;
   p3y << std::fixed << std::setprecision(3) << point3y;

   p4x << std::fixed << std::setprecision(3) << point4x;
   p4y << std::fixed << std::setprecision(3) << point4y;
   
   //String modification
   std::string line1 = "		" + p1x.str() + " " + p1y.str();
   std::string line2 = "		" + p2x.str() + " " + p2y.str();
   std::string line3 = "		" + p3x.str() + " " + p3y.str();
   std::string line4 = "		" + p4x.str() + " " + p4y.str();

   //Write down the 4 polygons and adjust location
    for(int i = 2; i < 6; i++) {
        double x = 1.95;
        double z = 1.95;
        if(i == 3 || i == 5) {
			x = x + 0.39;
		}
		if(i > 3 ) {
			z = z + 0.39;
		}
        outfile << "polygon " + std::to_string(i) << std::endl;
        outfile << "	count = 5" << std::endl;
        outfile << "	points =" << std::endl;

        outfile << line1 << std::endl;
        outfile << line2 << std::endl;
        outfile << line3 << std::endl;
        outfile << line4 << std::endl;
        outfile << line1 << std::endl;

        outfile << "	end points" << std::endl;
        outfile << "	angle = -180" << std::endl;

        std::stringstream sx;
        std::stringstream sz;

        sx << std::fixed << std::setprecision(2) << x;
        sz << std::fixed << std::setprecision(2) << z;

        outfile << "	begin.x = -5*Ax+-5*Cx+OffsetX+" + sx.str() << std::endl;
        outfile << "	begin.z = -5*Az+-5*Cz+OffsetZ+" + sz.str() << std::endl;
        outfile << "	begin.height = 2*Radius" << std::endl;
        outfile << "end polygon" << std::endl;
        outfile << "" << std::endl;
    }

    /* Read and write TextBlock.txt */
    if (textBlock.is_open()) {
        if (outfile.is_open()) {
            std::string line;
            while (std::getline(textBlock, line)) {
                outfile << line << std::endl;
            }
        }
        else {
            std::cerr << "Error: Could not open output file." << std::endl;
        }
    }
    else {
        std::cerr << "Unable to open textBlock.txt." << std::endl;
    }
    textBlock.close();
    outfile.close();    
}

/* ② */
void callPwe(){
    std::string command = bandSolvePath + " -c " + indPath;
    int result = system(command.c_str());
    if (result == 0) {
        std::cout << "Plane-wave Expansion executed successfully." << std::endl;
    } else {
        std::cout << "Plane-wave Expansion failed to execute." << std::endl;
    }
}

/* ③ */
void loadPweData(){
    std::ifstream pweData(pweDataPath);    
    std::string line;
    std::getline(pweData, line); //Skip the first line
    std::getline(pweData, line); //Read the second line

    std::string tokens[6]; //6 tokens
    std::string delimiter = "   ";
    size_t pos = 0;

    int i = 0;
    //First 5 tokens
    while ((pos = line.find(delimiter)) != std::string::npos && i < 5) {
        std::string token = line.substr(0, pos);
        tokens[i] = token;
        line.erase(0, pos + delimiter.length());
        i++;
    }
    //The last token
    tokens[5] = line.substr(0, 17);

    //Assign value to ωA,ωB,ωC,ωD
    try {        
        wA = std::stod(tokens[2]);
        wB = std::stod(tokens[3]);
        wC = std::stod(tokens[4]);
        wD = std::stod(tokens[5]);

        // std::cout << tokens[2] << std::endl;
        // std::cout << tokens[3] << std::endl;
        // std::cout << tokens[4] << std::endl;
        // std::cout << tokens[5] << std::endl;
        
    } catch (std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

/* ④ */
void processData() {
	double period = p1 * 1000;                 //Period(287.0)
	double a = period * pow(10, -7);           //a
	double b = (2 * M_PI) / a;                 //β0		
	double c = 2.99792458 * pow(10, 10);       //speed of light in cm/s
		
	double convertWA = wA * (2 * M_PI * c) / a;
	double convertWB = wB * (2 * M_PI * c) / a;
	double convertWC = wC * (2 * M_PI * c) / a;
	double convertWD = wD * (2 * M_PI * c) / a;

    double nav = (b*c)*(1+convertWA/convertWB)*2/(convertWA+convertWB+convertWC+convertWD);
    k3 = fabs(b - (wB * 2 * M_PI * c / ( a ))*nav/c);//fabs: float absolute value
    k1 = sqrt((b-((convertWC*nav)/c)+k3)*(b-k3)/4);
	k2d = ( 2 * pow(k1, 2))/b;
    k1d = k3;
    
    std::cout << "β0  : " << b  << std::endl;
    std::cout << "κ1  : " << k1 << std::endl;
    std::cout << "κ3  : " << k3 << std::endl;
    std::cout << "κ1D : " << k1d << std::endl;
    std::cout << "κ2D : " << k2d << std::endl;

    p2 = k1d;
    p3 = k2d;
}

/* ⑤ */
void callPmc(){
    std::string command =  exePath + " " 
    + std::to_string(p1) + " " + std::to_string(p2) + " " 
    + std::to_string(p3) + " " + std::to_string(p4) + " " 
    + std::to_string(p5) + " " + std::to_string(p6) + " "
    + std::to_string(p7) + " " + p8;
    
    int result = system(command.c_str());

    //Command execution test:
    if (result == 0) {
        std::cout << "Probabilistic Markov Chain executed successfully." << std::endl;
    } else {
        std::cout << "Probabilistic Markov Chain failed to execute." << std::endl;
    }  
}

/* ⑥ */
void recordAlphaParallel(){
    std::ofstream data(recordPath); //New output data recording file
    std::ifstream csv(p8);          //csv file generated from last step

    data << "No.,wA,wB,wC,wD,K1,K3,K1D,K2D,a||\n" + std::to_string(num) + ","; //Headers
    data << std::to_string(wA)+","+std::to_string(wB)+","+std::to_string(wC)+","+std::to_string(wD)+",";
    data << std::to_string(k1)+","+std::to_string(k3)+","+std::to_string(k1d)+","+std::to_string(k2d)+",";

    // i/o Alpha Parallel
    if (csv.is_open()) {     
        std::string line;
        std::getline(csv, line);
        std::string ap = line.substr(line.length() - 15);
        data << ap << std::endl;  
    }
    else {
        std::cout << "Unable to open csv file." << std::endl;
    }
    csv.close();
    data.close();
}

/* Main */
int main() {  
   createIndFile();         /* ① */
   callPwe();               /* ② */ 
   loadPweData();           /* ③ */
   processData();           /* ④ */
   callPmc();               /* ⑤ */
   recordAlphaParallel();   /* ⑥ */
   return 0;
}