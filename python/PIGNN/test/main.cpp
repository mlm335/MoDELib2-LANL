#include <pybind11/embed.h>
#include <filesystem>
#include <iostream>
#include <pybind11/eigen.h>
#include <string>
#include <pybind11/stl.h>
#include <list>
#include <map>
//#include <DDconfigIO.h>
//#include <AtomDisplacer.h>

namespace fs = std::filesystem;
namespace py = pybind11;

// Turn the GNN into a class, store trained model as constructor
class dglContainer
{
    public:
    py::object net;
//    py::module mtest = py::module::import("testHetetro");
    
    dglContainer()
    {
        std::cout << "Constructing the dgl" << std::endl;
        py::module mtest = py::module::import("MLmobility");
//        mtest = py::module::import("testHetetro");
        net = mtest.attr("PIGNN")();
    }
    
    void printforward(){
        std::cout << net.attr("read_testdata")() << std::endl;
    };
    
};


int main() {
    
    py::scoped_interpreter guard{};
    py::module sys = py::module::import("sys");
    py::list path = sys.attr("path");
    path.append("/Users/matthewmaron/Documents/Embedding/pybind7/");
    
    
    // Object that contains trained ML model:
    dglContainer myNet;
    myNet.printforward();
    //myNet.printforward();
    //myNet.printforward();
}
