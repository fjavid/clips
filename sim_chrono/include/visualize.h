#define PICOPARAVIEW_IMPLEMENTATION
#include "picoparaview.h"
// #include "include/system.h"

// #include <random>
// #include <chrono>
#include <Eigen/Core>

struct InfoVTU
{
    InfoVTU(double _time, int _ts, int _pnum, std::string _fn) : time(_time), ts_num(_ts), part_num(_pnum), file_name(_fn) {};
    double time;
    int ts_num;
    int part_num;
    std::string file_name;
};

class ParaVisual
{
    public:
    ParaVisual(std::string _path, std::string _name, double _dt) : path(_path), name(_name), dt(_dt) {};
    // ~ParaVisal(){};
    std::string path;
    std::string name;
    double dt;
    std::vector<InfoVTU> info_vtu; 

    void write_pvd(std::string pvd_fname)
    {
        // using namespace ppv::fmt;
        std::ofstream pvd(path + "/" + pvd_fname + ".pvd");
        ppv::fmt::OpenPVD( pvd );
        {
            for (auto ivtu : info_vtu)
                ppv::fmt::RecordPVDDataSet(pvd, ivtu.file_name, ivtu.time);
        }
        ppv::fmt::ClosePVD( pvd );
    }

    void write_wire_obj(std::vector<std::vector<double>> vertices, std::vector< std::vector<int>> edges, std::string obj_name, int id, int ts)
    {
        std::string f_name = obj_name + "_" + std::to_string(id) + "_t_" + std::to_string(ts) + ".vtu";
        double time = ts * dt;
        info_vtu.push_back(InfoVTU(time, ts, id, f_name));
        int nv = vertices.size();
        int ne = edges.size();
        ppv::fmt::e_VtkType elemType = ppv::fmt::LINE;
        std::ofstream vtu(path + '/' + f_name);
        ppv::fmt::OpenVTU(vtu, nv, ne);
        ppv::WriteGeometryVTU(vtu, vertices, edges, elemType);
        ppv::fmt::CloseVTU(vtu);
        vtu.close();
    }

    void write_face_obj(std::vector<std::vector<double>> vertices, std::vector< std::vector<int>> faces, std::string obj_name, int id, int ts)
    {
        std::string f_name = obj_name + "_" + std::to_string(id) + "_t_" + std::to_string(ts) + ".vtu";
        double time = ts * dt;
        info_vtu.push_back(InfoVTU(time, ts, id, f_name));
        int nv = vertices.size();
        int ne = faces.size();
        ppv::fmt::e_VtkType elemType = ppv::fmt::TRI;
        std::ofstream vtu(path + '/' + f_name);
        ppv::fmt::OpenVTU(vtu, nv, ne);
        ppv::WriteGeometryVTU(vtu, vertices, faces, elemType);
        ppv::fmt::CloseVTU(vtu);
        vtu.close();
    }

};