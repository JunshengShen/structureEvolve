#pragma once

//morton tool class
class suMorton
{
public:
	static int encode(int x, int y, int z, int level)
	{
		int morton = 0;
		for (int i = level - 1; i >= 0; i--)
		{
			morton += (x >> i) << (3 * i + 2);
			x = x - ((x >> i) << i);
			morton += (y >> i) << (3 * i + 1);
			y = y - ((y >> i) << i);
			morton += (z >> i) << (3 * i);
			z = z - ((z >> i) << i);
		}
		return morton;
	}
	static void decode(int &x, int &y, int &z, int morton, int level)
	{
		x = y = z = 0;
		for (int i = level - 1; i >= 0; i--)
		{
			x += (morton >> (i * 3 + 2)) << i;
			morton -= (morton >> (i * 3 + 2)) << (i * 3 + 2);
			y += (morton >> (i * 3 + 1)) << i;
			morton -= (morton >> (i * 3 + 1)) << (i * 3 + 1);
			z += (morton >> (i * 3)) << i;
			morton -= (morton >> (i * 3)) << (i * 3);
		}
		return;
	}

	static int get_6neighbors(std::vector<int> &neighbors, int morton_code, int level)  //return 6 neighbors
	{
		int nn[6][3] = {
			{ 1, 0, 0},
			{ -1, 0, 0},
			{ 0, 1, 0},
			{ 0, -1, 0},
			{ 0, 0, 1},
			{ 0, 0, -1}
		};
		int x, y, z;
		decode(x, y, z, morton_code, level);
		std::cout << x << ' ' << y << ' ' << z << std::endl;
		neighbors.clear();
		if (x == 0)
			nn[1][0] = 0;
		if (x == (1 << level)-1)
			nn[0][0] = 0;

		if (y == 0)
			nn[3][1] = 0;
		if (y == (1 << level)-1)
			nn[2][1] = 0;

		if (z == 0)
			nn[5][2] = 0;
		if (z == (1 <<  level)-1)
			nn[4][2] = 0;
		for (int i = 0; i < 6; i++)
		{
			neighbors.push_back(encode(x + nn[i][0], y + nn[i][1], z + nn[i][2], level) );
		}
		for (int i = 5; i >=0; i--) {
			if (neighbors[i] == morton_code)
				neighbors.erase(neighbors.begin()+i);
		}
		std::cout << neighbors.size();
		return neighbors.size();
	}
	static int get_26neighbors(std::vector<int> &neighbors, int morton_code, int level) {//return 26 neighbors
		int nn[26][3] = {
			{ 1, 0, 0 },
			{ -1, 0, 0 },
			{ 0, 1, 0 },
			{ 0, -1, 0 },
			{ 0, 0, 1 },
			{ 0, 0, -1 },
			{1,1,0},
			{-1,1,0},
			{-1,-1,0},
			{1,-1,0},
			{1,0,1},
			{-1,0,1},
			{1,0,-1},
			{-1,0,-1},
			{0,1,1},
			{0,-1,1},
			{0,-1,-1},
			{0,1,-1},
			{1,1,1},
			{-1,-1,-1},
			{1,1,-1},
			{1,-1,1},
			{-1,1,1},
			{1. - 1,-1},
			{-1,1,-1},
			{-1,-1,1}
		};
		int x, y, z;
		decode(x, y, z, morton_code, level);
		//std::cout << x << ' ' << y << ' ' << z << std::endl;
		neighbors.clear();

		//boundary test
		if (x == 0) {
			for (int i = 0; i < 26; i++) {
				if (nn[i][0] == -1) nn[i][0] = 0;
			}
		}

		if (x == (1 << level) - 1) {
			for (int i = 0; i < 26; i++) {
				if (nn[i][0] == 1) nn[i][0] = 0;
			}
		}

		if (y == 0) {
			for (int i = 0; i < 26; i++) {
				if (nn[i][1] == -1) nn[i][1] = 0;
			}
		}

		if (y == (1 << level) - 1) {
			for (int i = 0; i < 26; i++) {
				if (nn[i][1] == 1) nn[i][1] = 0;
			}
		}

		if (z == 0) {
			for (int i = 0; i < 26; i++) {
				if (nn[i][2] == -1) nn[i][2] = 0;
			}
		}

		if (z == (1 << level) - 1) {
			for (int i = 0; i < 26; i++) {
				if (nn[i][2] == 1) nn[i][2] = 0;
			}
		}


		for (int i = 0; i < 26; i++)
		{
			neighbors.push_back(encode(x + nn[i][0], y + nn[i][1], z + nn[i][2], level));
		}
		for (int i = 25; i >= 0; i--) {
			if (neighbors[i] == morton_code)
				neighbors.erase(neighbors.begin() + i);
		}
		std::cout << neighbors.size();
		return neighbors.size();

	}

};
