#include <igl/opengl/glfw/Viewer.h>
#include "io.h"
#include "gui.h"

using namespace Eigen;
using namespace std;
using namespace mvis;


int main(int argc, char *argv[])
{
    vector<MatrixXd> V_list;
    vector<MatrixXi> F_list;
    GUI gui = GUI(V_list, F_list);
    gui.load_example();
    gui.show();
}
