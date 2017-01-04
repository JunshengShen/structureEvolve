#include "overlap.h"
#include <vector>
int overlap(double boxcenter[3], double boxhalfsize[3], double triverts[3][3])
{
	double test1 = boxcenter[0], test2 = boxcenter[1], test3 = boxcenter[2];
	double test21 = boxhalfsize[0], test22 = boxhalfsize[1], test23 = boxhalfsize[2];
	double test31= triverts[0][0], test32 = triverts[0][1], test33 = triverts[0][2],
		test34 = triverts[1][0], test35 = triverts[1][1], test36 = triverts[1][2],
		test37 = triverts[2][0], test38 = triverts[2][1], test39 = triverts[2][2];
	double xMin = boxcenter[0] - boxhalfsize[0];
	double xMax = boxcenter[0] + boxhalfsize[0];
	double yMin = boxcenter[1] - boxhalfsize[1];
	double yMax = boxcenter[1] + boxhalfsize[1];
	double zMin = boxcenter[2] - boxhalfsize[2];
	double zMax = boxcenter[2] + boxhalfsize[2];
	std::vector<double*> insidePoints;
	std::vector<double*> outsidePoints;
	for (int i = 0; i < 3; i++) {//���м����ڵ����ı�����
		if (triverts[i][0] >= xMin&&triverts[i][0] <= xMax)
			if (triverts[i][1] >= yMin&&triverts[i][1] <= yMax)
				insidePoints.push_back(triverts[i]);
			else
				outsidePoints.push_back(triverts[i]);
		else
			outsidePoints.push_back(triverts[i]);
	}
	if(!insidePoints.size())//ȫ���ڵ����ı�����   ������
		return 0;
	for (int i = 0; i < insidePoints.size(); i++) {//����ж���������ͽ���
		if (insidePoints[i][2] >= zMin&&insidePoints[i][2] <= zMax)
			return 1;
	}
	if (insidePoints.size() == 3) {//ȫ��������
		int count = 0;
		for (int i = 0; i < 3; i++)
			if (insidePoints[i][2] > zMax)
				count++;
		if (count == 1 || count == 2)
			return 1;//����е��������� �е���������  �ͽ���
		else
			return 0;//��Ȼ�Ͳ�����
	}
	if (insidePoints.size() == 2) {//������������
		int count = 0;
		for (int i = 0; i < 2; i++)
			if (insidePoints[i][2] > zMax)
				count++;
		if (count ==1)//��������һ������һ������ �ͽ���
			return 1;
		for (int i = 0; i < 2; i++) {
			double y, z;
			if (outsidePoints[0][0] > xMax) {//���ұ�
				y = insidePoints[i][1] - (insidePoints[i][1] - outsidePoints[0][1])*
					(insidePoints[i][0] - xMax) / (insidePoints[i][0] - outsidePoints[0][0]);
				if (y > yMax || y < yMin) {//��X����������߽���
					if (outsidePoints[0][1 > yMax]) {//��ǰ��
						z = insidePoints[i][2] - (insidePoints[i][2] - outsidePoints[0][2])*
							(insidePoints[i][1] - yMax) / (insidePoints[i][1] - outsidePoints[0][1]);
						if (z > zMin&&z < zMin)
							return 1;
						else
							return 0;
					}

					else {//�ں���
						z = insidePoints[i][2] - (insidePoints[i][2] - outsidePoints[0][2])*
							(insidePoints[i][1] - yMin) / (insidePoints[i][1] - outsidePoints[0][1]);
						if (z > zMin&&z < zMin)
							return 1;
						else
							return 0;
					}
				}
				else {//��Y����������߽���
					z=insidePoints[i][2]- (insidePoints[i][2] - outsidePoints[0][2])*
						(insidePoints[i][0] - xMax) / (insidePoints[i][0] - outsidePoints[0][0]);
					if (z >= zMin&&z <= zMax)
						return 1;
					else return 0;
				}
			}
			else if (outsidePoints[0][0] < xMin) {//�����
				y = insidePoints[i][1] - (insidePoints[i][1] - outsidePoints[0][1])*
					(insidePoints[i][0] - xMin) / (insidePoints[i][0] - outsidePoints[0][0]);
				if (y > yMax || y < yMin) {//��X����������߽���
					if (outsidePoints[0][1 > yMax]) {//��ǰ��
						z = insidePoints[i][2] - (insidePoints[i][2] - outsidePoints[0][2])*
							(insidePoints[i][1] - yMax) / (insidePoints[i][1] - outsidePoints[0][1]);
						if (z > zMin&&z < zMin)
							return 1;
						else
							return 0;
					}

					else {//�ں���
						z = insidePoints[i][2] - (insidePoints[i][2] - outsidePoints[0][2])*
							(insidePoints[i][1] - yMin) / (insidePoints[i][1] - outsidePoints[0][1]);
						if (z > zMin&&z < zMin)
							return 1;
						else
							return 0;
					}
				}
				else {//��Y����������߽���
					z = insidePoints[i][2] - (insidePoints[i][2] - outsidePoints[0][2])*
						(insidePoints[i][0] - xMax) / (insidePoints[i][0] - outsidePoints[0][0]);
					if (z >= zMin&&z <= zMax)
						return 1;
					else return 0;
				}

			}
			
			else if (outsidePoints[0][1 > yMax]) {//��ǰ��
				z = insidePoints[i][2] - (insidePoints[i][2]-outsidePoints[0][2])*
					(insidePoints[i][1]-yMax) / (insidePoints[i][1]-outsidePoints[0][1]);
				if (z > zMin&&z < zMin)
					return 1;
				else
					return 0;
			}

			else {//�ں���
				z = insidePoints[i][2] - (insidePoints[i][2] - outsidePoints[0][2])*
					(insidePoints[i][1] - yMin) / (insidePoints[i][1] - outsidePoints[0][1]);
				if (z > zMin&&z < zMin)
					return 1;
				else
					return 0;
			}

		}



	}



	if (insidePoints.size() == 1) {
		for (int i = 0; i < 2; i++) {
			double y, z;
			if (outsidePoints[i][0] > xMax) {//���ұ�
				y = insidePoints[0][1] - (insidePoints[0][1] - outsidePoints[i][1])*
					(insidePoints[0][0] - xMax) / (insidePoints[0][0] - outsidePoints[i][0]);
				if (y > yMax || y < yMin) {//��X����������߽���
					if (outsidePoints[i][1 > yMax]) {//��ǰ��
						z = insidePoints[0][2] - (insidePoints[0][2] - outsidePoints[i][2])*
							(insidePoints[0][1] - yMax) / (insidePoints[0][1] - outsidePoints[i][1]);
						if (z > zMin&&z < zMin)
							return 1;
						else
							return 0;
					}

					else {//�ں���
						z = insidePoints[0][2] - (insidePoints[0][2] - outsidePoints[i][2])*
							(insidePoints[0][1] - yMin) / (insidePoints[0][1] - outsidePoints[i][1]);
						if (z > zMin&&z < zMin)
							return 1;
						else
							return 0;
					}
				}
				else {//��Y����������߽���
					z = insidePoints[0][2] - (insidePoints[0][2] - outsidePoints[i][2])*
						(insidePoints[0][0] - xMax) / (insidePoints[0][0] - outsidePoints[i][0]);
					if (z >= zMin&&z <= zMax)
						return 1;
					else return 0;
				}
			}
			else if (outsidePoints[i][0] < xMin) {//�����
				y = insidePoints[0][1] - (insidePoints[0][1] - outsidePoints[i][1])*
					(insidePoints[0][0] - xMin) / (insidePoints[0][0] - outsidePoints[i][0]);
				if (y > yMax || y < yMin) {//��X����������߽���
					if (outsidePoints[i][1 > yMax]) {//��ǰ��
						z = insidePoints[0][2] - (insidePoints[0][2] - outsidePoints[i][2])*
							(insidePoints[0][1] - yMax) / (insidePoints[0][1] - outsidePoints[i][1]);
						if (z > zMin&&z < zMin)
							return 1;
						else
							return 0;
					}

					else {//�ں���
						z = insidePoints[0][2] - (insidePoints[0][2] - outsidePoints[i][2])*
							(insidePoints[0][1] - yMin) / (insidePoints[0][1] - outsidePoints[i][1]);
						if (z > zMin&&z < zMin)
							return 1;
						else
							return 0;
					}
				}
				else {//��Y����������߽���
					z = insidePoints[0][2] - (insidePoints[0][2] - outsidePoints[i][2])*
						(insidePoints[0][0] - xMax) / (insidePoints[0][0] - outsidePoints[i][0]);
					if (z >= zMin&&z <= zMax)
						return 1;
					else return 0;
				}

			}

			else if (outsidePoints[i][1 > yMax]) {//��ǰ��
				z = insidePoints[0][2] - (insidePoints[0][2] - outsidePoints[i][2])*
					(insidePoints[0][1] - yMax) / (insidePoints[0][1] - outsidePoints[i][1]);
				if (z > zMin&&z < zMin)
					return 1;
				else
					return 0;
			}

			else {//�ں���
				z = insidePoints[0][2] - (insidePoints[0][2] - outsidePoints[i][2])*
					(insidePoints[0][1] - yMin) / (insidePoints[0][1] - outsidePoints[i][1]);
				if (z > zMin&&z < zMin)
					return 1;
				else
					return 0;
			}

		}
	}




}