#include "../src/GLE.h"
#include <vector>
#include <math.h>

int main(int argc, char* argv[])
{
     std::vector<double> x;
     std::vector<double> y1;
     std::vector<double> y2;
     std::vector< std::vector<double> > all_y;
     for (double i; i < 10; ++i)
    {
         x.push_back(i);
         y1.push_back(pow(i,2));
         y2.push_back(pow(i,3));
    }

    all_y.push_back(y1);
    all_y.push_back(y2);

    GLE g;
    GLE::PanelID panelID;
    GLE::PlotProperties plotProperties;
    GLE::PanelProperties panelProperties;
    panelProperties.title = "GLE++ Example";
    panelID = g.plot(x, all_y, plotProperties);
    g.setPanelProperties(panelProperties, panelID);
    g.canvasProperties.width = 6;
    g.canvasProperties.height = 6;
    g.draw("output.pdf");

    return 0;
}
