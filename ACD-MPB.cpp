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

/* Must activate mpb virtual environment before using */

int num = 1;               //Testing number

//Parameters for PMC
double wA, wB, wC, wD;     //Frequencies
double k1, k3, k1d, k2d;   //Kappas
double p1 = 0.287;         //period(mm) double
double p2 = k1d;           //K1D double (1334.3984)
double p3 = k2d;           //K2D double (386.8075)
double p4 = 38.0;          //Krad double
double p5 = 0.0;           //Ki double
double p6 = 0.0;           //reflectivity double
int p7 = 50;               //Time-step

//File storing paths:
//PLEASE ENTER YOUR LOCAL PATH HERE:
std::string exePath     = "/Users/.../ACD/bin/PMC.exe";    //Probablistic Markov Chain executable
std::string recordPath  = "/Users/.../ACD/data.csv";
std::string p8          = "/Users/.../ACD/PMC.csv";
std::string genePath    = "/Users/.../ACD/gene.txt";
std::string scriptPath  = "/Users/.../ACD/script.py";
std::string headerPath  = "/Users/.../ACD/header.txt"; //Scripts before creating cubes
std::string footerPath  = "/Users/.../ACD/footer.txt"; //Scripts after creating cubes
std::string pweDataPath = "/Users/.../ACD/script.dat"; //The auto-generated script


//PWE commands
std::string command1 = "python script.py > script.out";       //Run MPB simulation
std::string command2 = "grep freqs: script.out > script.dat"; //Grab useful data


//=== Dashboard ===
int n = 9;               //Grid resolution (must be an odd number)


int m = (n+1)/2;         //Median number
int arraySize  = n * n;


/* ① */
void geneGenerator(){
    /* Step 1 */
    // Seed the random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Define the distribution (1 or 0)
    std::uniform_int_distribution<int> dist(0, 1);

    // Generate an array of random 1s and 0s
    std::vector<int> geneArray(arraySize);
    for (int i = 0; i < arraySize; ++i) {
        geneArray[i] = dist(gen);
    }

    // Print the generated array to a txt file (as a record)
    std::ofstream txt(genePath);
    if (txt.is_open()) {
        for (int i : geneArray) {
            txt << i << " ";
        }
    }
    else {
            std::cerr << "Error: Could not open txt file." << std::endl;
    }

    /* Step 2 */
    // To write a python script for MPB   
    std::ofstream py(scriptPath);
    std::ifstream header(headerPath);
    std::ifstream footer(footerPath);

    /* Read and write header.txt line by line into outfile */
    if (header.is_open()) {
        if (py.is_open()) {
            std::string line;
            while (std::getline(header, line)) {
                py << line << std::endl;
            }
        }
        else {
            std::cerr << "Error: Could not open py file." << std::endl;
        }
    }
    else {
        std::cerr << "Unable to open header.txt." << std::endl;
    }
    header.close();

    //Parameters (Do calculations inside python)
    py << "n  = " + std::to_string(n) + "       # Grid resolution" << std::endl; //
    py << "a  = 287     # Period"          << std::endl;
    py << "ul = 175     # Unit Length"     << std::endl;
    py << "gl = ul/a    # Grid Length"     << std::endl;
    py << "pl = gl/n    # Pixel Length"    << std::endl;
    py << "v  = 1/2*pl  # Vector"          << std::endl;
    py << "vertices = [mp.Vector3(v,v),mp.Vector3(-v,v),mp.Vector3(-v,-v),mp.Vector3(v,-v)]" << std::endl;
    py << "geometry = [";

    std::string prism = "mp.Prism(vertices, height=mp.inf, center=mp.Vector3(";
    std::string material = "),material=mp.Medium(epsilon=1)),";  //Cubes are airholes (epsilon=1)

    //Traverse the array and process
    for (int i = 0; i < arraySize; ++i) {
        if (geneArray[i] == 0) {      //If the Gene fragment is 0, skip
            continue; // Skip
        }
        else if (geneArray[i] == 1){  //If the Gene fragment is 1, create the cube, and place it
            int pn = i + 1;           //position number (arrays start from 0, position is therefore +1)
            int rowNum = ceil(pn/n);  //Calculate its row number
            int colNum = pn%n;        //Calculate its column number

            std::string rowStr = std::to_string(rowNum-m); // y displacement
            std::string colStr = std::to_string(colNum-m); // x displacement

            py << prism + colStr + "*pl," + rowStr + "*pl" + material;
        }       
    }
    py << "]" << std::endl;

    /* add footer.txt */
    if (footer.is_open()) {
        if (py.is_open()) {
            std::string line;
            while (std::getline(footer, line)) {
                py << line << std::endl;
            }
        }
        else {
            std::cerr << "Error: Could not open py file." << std::endl;
        }
    }
    else {
        std::cerr << "Unable to open footer.txt." << std::endl;
    }
    footer.close();
    py.close();    
}

/* ② */
void callPwe(){    
    int result1 = system(command1.c_str());            // Command1 Run MPB
    if (result1 == 0) {
        std::cout << "Plane-wave Expansion executed successfully." << std::endl;
    } else {
        std::cout << "Plane-wave Expansion failed to execute." << std::endl;
    }

    int result2 = system(command2.c_str());            // Command2 Grab data
}

/* ③ */
void loadPweData(){

    std::ifstream pweData(pweDataPath);    
    std::string line;
    std::getline(pweData, line);  //Skip the 1st line
    std::getline(pweData, line);  //Skip the 2nd line
    std::getline(pweData, line);  //Read the 3rd line

    std::string tokens[11];        //11 tokens
    std::string delimiter = ", ";
    size_t pos = 0;

    int i = 0;
    //First 10 tokens
    while ((pos = line.find(delimiter)) != std::string::npos && i < 10) {
        std::string token = line.substr(0, pos);
        tokens[i] = token;
        line.erase(0, pos + delimiter.length());
        i++;
    }
    //The last token
    tokens[10] = line.substr(0, 8);

    //Assign values to ωA,ωB,ωC,ωD
    try {
        wA = std::stod(tokens[7]);
        wB = std::stod(tokens[8]);
        wC = std::stod(tokens[9]);
        wD = std::stod(tokens[10]);

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

    if (result == 0) {
        std::cout << "Probabilistic Markov Chain executed successfully." << std::endl;
    } else {
        std::cout << "Probabilistic Markov Chain failed." << std::endl;
    }  
}

/* ⑥ */
void recordAlphaParallel(){
    std::ofstream data(recordPath); //New data-recording output file
    std::ifstream csv(p8);          //csv file generated from Markov Chain

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
   geneGenerator();         /* ① */
   callPwe();               /* ② */ 
   loadPweData();           /* ③ */
   processData();           /* ④ */
   callPmc();               /* ⑤ */
   recordAlphaParallel();   /* ⑥ */
   return 0;
}