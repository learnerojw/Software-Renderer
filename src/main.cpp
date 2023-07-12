#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
//#include<algorithm>
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
Model* model = NULL;
const int width = 800;
const int height = 800;
Vec3f light_dir(0, 0, -1);
Vec3f camera(0, 0, 3);
const int depth = 255;

//叉乘
Vec3f Cross(Vec3f v1, Vec3f v2)
{
    return Vec3f(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}
void line(Vec2i p0, Vec2i p1, TGAImage& image, TGAColor color)
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

    for (int x = p0.x; x <= p1.x; x++)
    {
        float t = (x - p0.x) / float(p1.x - p0.x);
        int y = p0.y * (1 - t) + p1.y * t;
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

//判断屏幕空间的像素点是否在屏幕空间的一个三角形内
bool insideTriangle(Vec3f t0, Vec3f t1, Vec3f t2, Vec3f p)
{
    Vec3f t0_t1 = t1 - t0;
    Vec3f t0_p = p - t0;

    Vec3f t1_t2 = t2 - t1;
    Vec3f t1_p = p - t1;

    Vec3f t2_t0 = t0 - t2;
    Vec3f t2_p = p - t2;

    Vec3f res1 = Vec3f(t0_t1.y * t0_p.z - t0_p.y * t0_t1.z, t0_t1.z * t0_p.x - t0_p.z * t0_t1.x, t0_t1.x * t0_p.y - t0_p.x * t0_t1.y);
    Vec3f res2 = Vec3f(t1_t2.y * t1_p.z - t1_p.y * t1_t2.z, t1_t2.z * t1_p.x - t1_p.z * t1_t2.x, t1_t2.x * t1_p.y - t1_p.x * t1_t2.y);
    Vec3f res3 = Vec3f(t2_t0.y * t2_p.z - t2_p.y * t2_t0.z, t2_t0.z * t2_p.x - t2_p.z * t2_t0.x, t2_t0.x * t2_p.y - t2_p.x * t2_t0.y);

    return (res1.z >= 0 && res2.z >= 0 && res3.z >= 0) || (res1.z <= 0 && res2.z <= 0 && res3.z <= 0);
}

//求重心坐标
Vec3f barycentric(Vec3f t0, Vec3f t1, Vec3f t2, Vec3f& p)
{
    t0.z = 0;
    t1.z = 0;
    t2.z = 0;
    Vec3f t0_t1 = t1 - t0;
    Vec3f t1_t2 = t2 - t1;
    Vec3f t2_t0 = t0 - t2;
    Vec3f p_t0 = t0 - p;
    Vec3f p_t1 = t1 - p;
    Vec3f p_t2 = t2 - p;

    Vec3f all = Vec3f(t0_t1.y * t1_t2.z - t1_t2.y * t0_t1.z, t1_t2.x * t0_t1.z - t1_t2.z * t0_t1.x, t0_t1.x * t1_t2.y - t1_t2.x * t0_t1.y);
    float all_size = sqrt(all.x * all.x + all.y * all.y + all.z * all.z);

    Vec3f res1 = Cross(p_t1, p_t2);
    float alpha = sqrt(res1.x * res1.x + res1.y * res1.y + res1.z * res1.z) / all_size;

    Vec3f res2 = Cross(p_t0, p_t2);
    float beta = sqrt(res2.x * res2.x + res2.y * res2.y + res2.z * res2.z) / all_size;

    Vec3f res3 = Cross(p_t0, p_t1);
    float gama = sqrt(res3.x * res3.x + res3.y * res3.y + res3.z * res3.z) / all_size;

    return Vec3f(alpha, beta, gama);
}


//主流方法绘制三角形
void triangle2(Vec3f t0, Vec3f t1, Vec3f t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, float* zbuffer, TGAImage& image, float intensity)
{
    //缩小遍历范围
    int x_min = floor(std::min(std::min(t0.x, t1.x), t2.x));
    int x_max = ceil(std::max(std::max(t0.x, t1.x), t2.x));

    int y_min = floor(std::min(std::min(t0.y, t1.y), t2.y));
    int y_max = ceil(std::max(std::max(t0.y, t1.y), t2.y));

    for (int i = x_min; i <= x_max; i++)
    {
        for (int j = y_min; j <= y_max; j++)
        {
            Vec3f p(i, j, 0);
            //判断该屏幕像素坐标是否在三角形内
            if (insideTriangle(t0, t1, t2, p))
            {
                //计算重心坐标
                Vec3f barycenter = barycentric(t0, t1, t2, p);
                //通过重心坐标对z值进行插值，这里还需优化
                float z = t0.z * barycenter.x + t1.z * barycenter.y + t2.z * barycenter.z;

                //通过重心坐标对uv坐标进行插值
                Vec2i uv = uv0 * barycenter.x + uv1 * barycenter.y + uv2 * barycenter.z;

                //深度测试
                if (z > zbuffer[j * width + i])
                {
                    //深度写入
                    zbuffer[j * width + i] = z;
                    //采样
                    TGAColor color = model->diffuse(uv);
                    image.set(i, j, TGAColor(color.r * intensity, color.g * intensity, color.b * intensity));
                }

            }
        }
    }
}

//获得三维坐标
Vec3f threePos(Matrix m)
{
    return Vec3f(m[0][0], m[1][0], m[2][0]);
}

//齐次坐标转换为普通坐标
Matrix m2v(Matrix m) {
    Matrix res(4, 1);
    res[0][0] = m[0][0] / m[3][0];
    res[1][0] = m[1][0] / m[3][0];
    res[2][0] = m[2][0] / m[3][0];
    res[3][0] = 1.0f;
    return res;
}

//普通三维坐标转换成齐次坐标
Matrix v2m(Vec3f v) {
    Matrix m(4, 1);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

//视口变换矩阵，这里跟入门精要里面的写法有区别
Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity(4);
    m[0][3] = x + w / 2.f;
    m[1][3] = y + h / 2.f;
    m[2][3] = depth / 2.f;

    m[0][0] = w / 2.f;
    m[1][1] = h / 2.f;
    m[2][2] = depth / 2.f;
    return m;
}


int main(int argc, char** argv) {
    if (2 == argc) {
        model = new Model(argv[1]);
    }
    else {
        model = new Model("obj/african_head.obj");
    }

    TGAImage image(width, height, TGAImage::RGB);

    //设置深度缓冲区大小
    float* zbuffer = new float[width * height * 2];
    for (int i = 0; i < width * height * 2; i++) {
        zbuffer[i] = std::numeric_limits<int>::min();
    }

    {
        Matrix Projection = Matrix::identity(4);
        Matrix ViewPort = viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
        //这里的投影平面设置为z=0的平面，可推导出该位置的值
        Projection[3][2] = -1.f / camera.z;

        //遍历每个面
        for (int i = 0; i < model->nfaces(); i++) {

            std::vector<int> face = model->face(i);
            Vec3i screen_coords[3];
            Vec3f world_coords[3];
            //遍历每个面的顶点
            for (int j = 0; j < 3; j++) {
                Vec3f v = model->vert(face[j]);
                //先转换为齐次坐标，然后进行（模型）（视图）投影矩阵的转换，然后再进行齐次除法得到NDC坐标，然后进行视口变换
                screen_coords[j] = threePos(ViewPort * m2v(Projection * v2m(v)));
                world_coords[j] = v;
            }
            //通过对三角形两条边进行叉乘，得到法线向量
            Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
            n.normalize();

            //光线与法线点乘
            float intensity = n * light_dir;
            if (intensity > 0) {
                Vec2i uv[3];
                for (int k = 0; k < 3; k++) {
                    uv[k] = model->uv(i, k);
                }
                triangle2(screen_coords[0], screen_coords[1], screen_coords[2], uv[0], uv[1], uv[2], zbuffer, image, intensity);
            }
        }

        image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
        image.write_tga_file("output.tga");
    }

    delete model;
    return 0;
}

