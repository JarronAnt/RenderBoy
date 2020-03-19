#pragma once
#include "vec3.h"
#include "aabb.h"


//this class repersents a node in our bounded volume heiarchy  
class bvhNode : public hittable {
	bvhNode() {}
	bvhNode(hittable **l, int n, float time0, float time1);

	virtual bool hit(const ray &r, float tmin, float tmax, hitRecord &hr) const;
	virtual bool boundingBox(float t0, float t1, aabb &box) const;

	hittable *left;
	hittable *right;
	aabb box;
};

bool bvhNode::boundingBox(float t0, float t1, aabb &b) const {
	b = box;
	return true;
}

//this function occurs when a node of our bvh gets hit
bool bvhNode::hit(const ray &r, float tmin, float tmax, hitRecord &hr) const {

	if (box.hit(r, tmin, tmax)) {
		hitRecord leftRec, rightRec;
		bool hitLeft = left->hit(r, tmin, tmax, leftRec);
		bool hitRight = right->hit(r, tmin, tmax, rightRec);

		//if we hit both the left and right store the closer one
		if (hitLeft && hitRight) {
			if (leftRec.t < rightRec.t) {
				hr = leftRec;
			}
			else {
				hr = rightRec;
			}
			return true;
		}
		//only hit left
		else if (hitLeft) {
			hr = leftRec;
			return true;
		}
		//only hit right
		else if (hitRight) {
			hr = rightRec;
		}
		else {
			return false;
		}

	}
	else 
		return false;
}


//this compares the boxes along the x axis
int box_x_compare(const void * a, const void * b) {
	aabb box_left, box_right;
	hittable *ah = *(hittable**)a;
	hittable *bh = *(hittable**)b;
	if (!ah->boundingBox(0, 0, box_left) || !bh->boundingBox(0, 0, box_right))
		std::cerr << "no bounding box in bvh_node constructor\n";
	if (box_left.min().x() - box_right.min().x() < 0.0)
		return -1;
	else
		return 1;
}

//this compares the boxes along the y axis
int box_y_compare(const void * a, const void * b)
{
	aabb box_left, box_right;
	hittable *ah = *(hittable**)a;
	hittable *bh = *(hittable**)b;
	if (!ah->boundingBox(0, 0, box_left) || !bh->boundingBox(0, 0, box_right))
		std::cerr << "no bounding box in bvh_node constructor\n";
	if (box_left.min().y() - box_right.min().y() < 0.0)
		return -1;
	else
		return 1;
}

//this compares the boxes along the z axis
int box_z_compare(const void * a, const void * b)
{
	aabb box_left, box_right;
	hittable *ah = *(hittable**)a;
	hittable *bh = *(hittable**)b;
	if (!ah->boundingBox(0, 0, box_left) || !bh->boundingBox(0, 0, box_right))
		std::cerr << "no bounding box in bvh_node constructor\n";
	if (box_left.min().z() - box_right.min().z() < 0.0)
		return -1;
	else
		return 1;
}

//this constructor builds our bvh
bvhNode::bvhNode(hittable **l, int n, float time0, float time1) {
	//get a random axis
	int axis = int(3 *  (rand() / (RAND_MAX + 1.0)));

	//run a quick sort on the list and compare based on the axis that was picked
	//note we use C qsort instead of C++ sort because we can explicitly set the comparison function 
	if (axis == 0) {
		qsort(l, n, sizeof(hittable*), box_x_compare);
	}
	else if (axis == 1) {
		qsort(l, n, sizeof(hittable*), box_y_compare);
	}
	else {
		qsort(l, n, sizeof(hittable*), box_z_compare);
	}

	//this part splits the list recusrively 
	if (n == 1) {
		left = right = l[0];
	}
	else if (n == 2) {
		left = l[0];
		right = l[1];
	}
	else {
		left = new bvhNode(l, n / 2, time0, time1);
		right = new bvhNode(l + n / 2, n - n / 2, time0, time1);
	}

	aabb leftBox, rightBox;

	if (!left->boundingBox(time0, time1, leftBox) || !right->boundingBox(time0, time1, rightBox)) {
		std::cerr << "no bounding box in bvh node constructor\n";
	}

	box = surroundingBox(leftBox, rightBox);
}