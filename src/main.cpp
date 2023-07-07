#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0, 255, 0, 255);
Model *model = NULL;
const int width  = 200;
const int height = 200;

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

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color)
{
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t0.y > t2.y) std::swap(t0, t2);
    if (t1.y > t2.y) std::swap(t1, t2);

    int total_height = t2.y - t0.y;
    int t0_t1_height = t1.y - t0.y + 1;
    int t1_t2_height = t2.y - t1.y + 1;
    for (int y = t0.y; y <= t2.y; y++)
    {
        float ooooo;
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


int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    TGAImage image(width, height, TGAImage::RGB);
    
    Vec2i t0[3] = { Vec2i(10,70),Vec2i(50,160),Vec2i(70,80) };
    Vec2i t1[3] = { Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180) };
    Vec2i t2[3] = { Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180) };
   
    triangle(t0[0], t0[1], t0[2], image, red);
    triangle(t1[0], t1[1], t1[2], image, white);
    triangle(t2[0], t2[1], t2[2], image, green);

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}

