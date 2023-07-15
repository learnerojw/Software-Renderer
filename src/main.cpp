#include <vector>
#include <iostream>

#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h"
//#include<algorithm>
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
Model* model = NULL;
const int width = 800;
const int height = 800;
Vec3f light_dir(1, 1, 1);
Vec3f       eye(1, 1, 3);
Vec3f    center(0, 0, 0);
Vec3f        up(0, 1, 0);
//���
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
//ɨ��ķ�������������
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

//�ж���Ļ�ռ�����ص��Ƿ�����Ļ�ռ��һ����������
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

//����������
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


//������������������
void triangle2(Vec2f* index, IShader& shader, TGAImage& zbuffer, TGAImage& image)
{
    Vec3f vertex[3];
    Vec4f screen_pos_ori[3];
    Vec3f screen_pos[3];
    TGAColor color;
    for (int i = 0; i < 3; i++)
    {
        /*std::vector<int> face = model->face(index[i][0]);
        vertex[i] = model->vert(face[index[i][1]]);*/
        screen_pos_ori[i]=shader.vertex(index[i].x, index[i].y);
        screen_pos[i] = Vec3f(screen_pos_ori[i][0] / screen_pos_ori[i][3], screen_pos_ori[i][1] / screen_pos_ori[i][3], screen_pos_ori[i][2] / screen_pos_ori[i][3]);
    }
    
    //��С������Χ
    int x_min = floor(std::min(std::min(screen_pos[0][0], screen_pos[1][0]), screen_pos[2][0]));
    //x_min = std::max(0, x_min);
    int x_max = ceil(std::max(std::max(screen_pos[0][0], screen_pos[1][0]),screen_pos[2][0]));
    //x_max = std::min(width, x_min);

    int y_min = floor(std::min(std::min(screen_pos[0][1], screen_pos[1][1]), screen_pos[2][1]));
    //y_min = std::max(0, y_min);
    int y_max = ceil(std::max(std::max(screen_pos[0][1], screen_pos[1][1]), screen_pos[2][1]));
    //y_max = std::min(height, y_max);

    for (int i = x_min; i <= x_max; i++)
    {
        for (int j = y_min; j <= y_max; j++)
        {
            Vec3f p(i, j, 0);
            //�жϸ���Ļ���������Ƿ�����������
            if (insideTriangle(screen_pos[0], screen_pos[1], screen_pos[2], p))
            {
                //������������
                Vec3f barycenter = barycentric(screen_pos[0], screen_pos[1], screen_pos[2], p);
                //ͨ�����������zֵ���в�ֵ�����ﻹ���Ż�
                float z = screen_pos[0].z * barycenter.x + screen_pos[1].z * barycenter.y + screen_pos[2].z * barycenter.z;

                //ͨ�����������uv������в�ֵ
                //Vec2i uv = uv0 * barycenter.x + uv1 * barycenter.y + uv2 * barycenter.z;
                int frag_depth = std::max(0, std::min(255, int(z+.5)));
                //Vec3f normal = normals[0] * barycenter.x + normals[1] * barycenter.y + normals[2] * barycenter.z;

                //float intensity = normal.normalize() * light_dir;
                //��Ȳ���
                if (zbuffer.get(i, j)[0] < frag_depth)
                {
                    //���д��
                    zbuffer.set(i, j, TGAColor(frag_depth));
                    shader.fragment(barycenter, color);
                    image.set(i, j, color);
                }

            }
        }
    }
}

//�����ά����
Vec3f threePos(Matrix m)
{
    return Vec3f(m[0][0], m[1][0], m[2][0]);
}

//�������ת��Ϊ��ͨ����
//Matrix m2v(Matrix m) {
//    Matrix res(4, 1);
//    res[0][0] = m[0][0] / m[3][0];
//    res[1][0] = m[1][0] / m[3][0];
//    res[2][0] = m[2][0] / m[3][0];
//    res[3][0] = 1.0f;
//    return res;
//}
//
////��ͨ��ά����ת�����������
//Matrix v2m(Vec3f v) {
//    Matrix m(4, 1);
//    m[0][0] = v.x;
//    m[1][0] = v.y;
//    m[2][0] = v.z;
//    m[3][0] = 1.f;
//    return m;
//}

//�ӿڱ任������������ž�Ҫ�����д��������
//Matrix viewport(int x, int y, int w, int h) {
//    Matrix m = Matrix::identity(4);
//    m[0][3] = x+w/2.f;
//    m[1][3] = y+h/2.0f;
//    m[2][3] = depth / 2.f;
//
//    m[0][0] = w / 2.f;
//    m[1][1] = h / 2.f;
//    m[2][2] = depth / 2.f;
//    return m;
//}

//ģ�͵��ӽǵ�ת��������û�п�������ռ�
//Matrix ModelView(Vec3f eye,Vec3f center,Vec3f up)
//{
//    Vec3f z = eye.normalize();
//    Vec3f x = Cross(up, eye).normalize();
//    Vec3f y = Cross(z, x).normalize();
//    Matrix res = Matrix::identity(4);
//
//    for (int i = 0; i < 3; i++)
//    {
//        res[0][i] = x[i];
//        res[1][i] = y[i];
//        res[2][i] = z[i];
//        res[i][3] = -center[i];
//    }
//
//    return res;
//}

struct GouraudShader : public IShader {
    Vec2f uvs[3];
    Vec3f normals[3];
    Matrix MVP = Projection * ModelView;
    Matrix MVP_invert = MVP.invert_transpose();
    virtual Vec4f vertex(int iface, int nthvert) {
        uvs[nthvert] = model->uv(iface, nthvert);
        normals[nthvert] = model->normal(iface, nthvert);
        Vec4f vertex = embed<4>(model->vert(iface, nthvert));
        Vec4f screen_pos= Viewport * Projection * ModelView * vertex;
        return screen_pos;
    }

    virtual bool fragment(Vec3f bar, TGAColor& color) {
        //Vec3f normal = normals[0] * bar.x + normals[1] * bar.y + normals[2] * bar.z;
        //normal.normalize();
        Vec2f uv = uvs[0] * bar.x + uvs[1] * bar.y + uvs[2] * bar.z;
        Vec3f normal= proj<3>(MVP_invert * embed<4>(model->normal(uv))).normalize();
        Vec3f l = proj<3>(MVP * embed<4>(light_dir)).normalize();
        Vec3f r = (normal * (normal * l * 2.f) - l).normalize();
        Vec3f viewDir= proj<3>(MVP * embed<4>(eye)).normalize();
        float diffuse_intensity = std::max(0.f, normal * l);
        float specular_intensity = pow(std::max(viewDir * r, 0.0f), 128);
        TGAColor diffuse = model->diffuse(uv) * diffuse_intensity;
        TGAColor specular = model->diffuse(uv) * specular_intensity*0.6f;
        
        color = diffuse + specular;
        return false;
    }
};

int main(int argc, char** argv) {
    if (2 == argc) {
        model = new Model(argv[1]);
    }
    else {
        model = new Model("obj/african_head.obj");
    }

    lookat(eye, center, up);
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    projection(-1.f / (eye - center).norm());
    light_dir.normalize();

    TGAImage image(width, height, TGAImage::RGB);
    TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

    GouraudShader shader;
    for (int i = 0; i < model->nfaces(); i++) {
        Vec2f index[3];
        for (int j = 0; j < 3; j++) {
            index[j]=Vec2f(i,j);
        }
        triangle2(index, shader, zbuffer, image);
    }

    image.flip_vertically(); // to place the origin in the bottom left corner of the image
    zbuffer.flip_vertically();
    image.write_tga_file("output.tga");
    zbuffer.write_tga_file("zbuffer.tga");

    delete model;
    return 0;
}

