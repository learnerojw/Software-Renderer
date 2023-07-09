#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
//#include<algorithm>
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0, 255, 0, 255);
Model *model = NULL;
const int width  = 800;
const int height = 800;

void line(Vec2i p0,Vec2i p1,TGAImage& image,TGAColor color)
{
    bool steep = false;

    if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y))
    {
        std::swap(p0.x, p0.y);
        std::swap(p1.x, p1.y);
        steep = true;
    }

    if (p0.x > p1.x)
    {
        std::swap(p0.x, p1.x);
        std::swap(p0.y, p1.y);
    }

    for (int x=p0.x;x<=p1.x;x++)
    {
        float t = (x - p0.x) / float(p1.x - p0.x);
        int y = p0.y * (1-t) + p1.y * t;
        if (steep)
        {
            image.set(y, x, color);
        }
        else
        {
            image.set(x, y, color);
        }
        
    }
}
//扫描的方法绘制三角形
void triangle1(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color)
{
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t0.y > t2.y) std::swap(t0, t2);
    if (t1.y > t2.y) std::swap(t1, t2);

    int total_height = t2.y - t0.y;
    int t0_t1_height = t1.y - t0.y + 1;
    int t1_t2_height = t2.y - t1.y + 1;
    for (int y = t0.y; y <= t2.y; y++)
    {
        float www;
        float alpha = (float)(y - t0.y) / total_height;
        float beta = (float)(y - t0.y) / t0_t1_height;
        float gama = (float)(y - t1.y) / t1_t2_height;


        Vec2i A = t0 + (t2 - t0) * alpha;
        Vec2i B;
        if (y <= t1.y)
        {
            B = t0 + (t1 - t0) * beta;
        }
        else
        {
            B = t1 + (t2 - t1) * gama;
        }

        image.set(A.x, y, color);
        image.set(B.x, y, color);

        if (A.x > B.x) std::swap(A, B);
        for (int x = A.x; x <= B.x; x++)
        {
            image.set(x, y, color);
        }
    }
}


bool insideTriangle(Vec3i t0, Vec3i t1, Vec3i t2,Vec3i p)
{
    Vec3i t0_t1 = t1 - t0;
    Vec3i t0_p = p - t0;

    Vec3i t1_t2 = t2 - t1;
    Vec3i t1_p = p - t1;

    Vec3i t2_t0 = t0 - t2;
    Vec3i t2_p = p - t2;
    
    Vec3i res1 = Vec3i(t0_t1.y * t0_p.z - t0_p.y * t0_t1.z, t0_t1.z * t0_p.x - t0_p.z * t0_t1.x, t0_t1.x * t0_p.y - t0_p.x * t0_t1.y);
    Vec3i res2 = Vec3i(t1_t2.y * t1_p.z - t1_p.y * t1_t2.z, t1_t2.z * t1_p.x - t1_p.z * t1_t2.x, t1_t2.x * t1_p.y - t1_p.x * t1_t2.y);
    Vec3i res3 = Vec3i(t2_t0.y * t2_p.z - t2_p.y * t2_t0.z, t2_t0.z * t2_p.x - t2_p.z * t2_t0.x, t2_t0.x * t2_p.y - t2_p.x * t2_t0.y);

    return (res1.z >= 0 && res2.z >= 0 && res3.z >= 0) || (res1.z <= 0 && res2.z <= 0 && res3.z <= 0);
}
//主流方法绘制三角形
void triangle2(Vec3i t0,Vec3i t1,Vec3i t2,TGAImage& image,TGAColor color)
{

    int x_min = std::min(std::min(t0.x, t1.x), t2.x);
    int x_max = std::max(std::max(t0.x, t1.x), t2.x);

    int y_min = std::min(std::min(t0.y, t1.y), t2.y);
    int y_max = std::max(std::max(t0.y, t1.y), t2.y);

    for (int i = x_min; i <= x_max; i++)
    {
        for (int j = y_min; j <= y_max; j++)
        {
            Vec3i p(i, j, 0);
            if (insideTriangle(t0, t1, t2, p))
            {
                image.set(i, j, color);
            }
        }
    }
}




int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    TGAImage image(width, height, TGAImage::RGB);
    
    /*Vec3i t0[3] = { Vec3i(10,70,0),Vec3i(50,160,0),Vec3i(70,80,0) };
    Vec3i t1[3] = { Vec3i(180, 50,0),  Vec3i(150, 1,0),   Vec3i(70, 180,0) };
    Vec3i t2[3] = { Vec3i(180, 150,0), Vec3i(120, 160,0), Vec3i(130, 180,0) };
   
    triangle2(t0[0], t0[1], t0[2], image, red);
    triangle2(t1[0], t1[1], t1[2], image, white);
    triangle2(t2[0], t2[1], t2[2], image, green);*/

    Vec3f light_dir(0, 0, -1);
    for (int i = 0; i < model->nfaces(); i++)//遍历所有面
    {
        std::vector<int> face = model->face(i);
        Vec3i screen_Pos[3];
        Vec3f obj_Poss[3];
        for (int j = 0; j < 3; j++)//遍历该面的所有顶点
        {
            Vec3f obj_Pos = model->vert(face[j]);
            screen_Pos[j] = Vec3i((obj_Pos.x + 1.0f) / 2.0f * width, (obj_Pos.y + 1.0f) / 2.0f * height, obj_Pos.z);//转换到屏幕空间
            obj_Poss[j] = obj_Pos;
        }
        Vec3f normal = (obj_Poss[2] - obj_Poss[0]) ^ (obj_Poss[1] - obj_Poss[0]);
        normal.normalize();
        float intensity = normal * light_dir;
        if (intensity > 0)
        {
            triangle2(screen_Pos[0], screen_Pos[1], screen_Pos[2], image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
        }
        
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}

