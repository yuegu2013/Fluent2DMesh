//
//  Face.h
//
//
//  Created by Ling Zou on 9/7/13.
//
//

#ifndef _Face_h
#define _Face_h

#include "Node.h"
#include "Cell.h"

class Face {
public:
  Face() : _node_id1(-1), _node_id2(-1), _cell_id1(-1), _cell_id2(-1),
  _id(-1) {}
  ~Face() {}

  Face(long int node_id1, long int node_id2, long int cell_id1, long int cell_id2, long int face_id) :
    _node_id1(node_id1), _node_id2(node_id2), _cell_id1(cell_id1), _cell_id2(cell_id2),
    _id(face_id)
  {}
  const long int node_id1() const { return _node_id1; }
  const long int node_id2() const { return _node_id2; }
  const long int cell_id1() const { return _cell_id1; }
  const long int cell_id2() const { return _cell_id2; }
  const long int id() const { return _id; }
  const double area() const { return _area; }
  const Vec3d & faceNormal() const { return _face_normal; }
  inline double & area() { return _area; }
  void reorderCells()
  {
    std::swap(_cell_id1, _cell_id2);
    // flip face normal
    //setFaceNormal(-_face_normal.x(), -_face_normal.y(), 0.0);
  }
  void setFaceNormal(double x, double y, double z)
  { _face_normal.x() = x; _face_normal.y() = y; _face_normal.z() = z; }

protected:
  long int _id, _node_id1, _node_id2, _cell_id1, _cell_id2;
  double _area;
  Vec3d _face_normal;
};

#endif
