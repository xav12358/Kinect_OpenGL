#include "eigen_func.h"

//////////////////////////////////////////////
/// \brief crossE
/// \param u
/// \param v
/// \return
///
Eigen::Vector3f crossE(Eigen::Vector3f u, Eigen::Vector3f v)
{
    Eigen::Vector3f result;

//    std::cout << "in u " << u << std::endl;
//    std::cout << "in v " << v << std::endl;

    result(0,0) = u(1,0)*v(2,0) - u(2,0)*v(1,0);
    result(1,0) = u(2,0)*v(0,0) - u(0,0)*v(2,0);
    result(2,0) = u(0,0)*v(1,0) - u(1,0)*v(0,0);

//    std::cout << "in result " << result << std::endl;
    return result;
}


////////////////////////////////////
/// \brief create_rotation_matrix
/// \param ax
/// \param ay
/// \param az
/// \return
///
Eigen::Affine3d create_rotation_matrix(double ax, double ay, double az) {
    Eigen::Affine3d rx =
            Eigen::Affine3d(Eigen::AngleAxisd(ax, Eigen::Vector3d(1, 0, 0)));
    Eigen::Affine3d ry =
            Eigen::Affine3d(Eigen::AngleAxisd(ay, Eigen::Vector3d(0, 1, 0)));
    Eigen::Affine3d rz =
            Eigen::Affine3d(Eigen::AngleAxisd(az, Eigen::Vector3d(0, 0, 1)));
    return rx*ry*rz;
}


/////////////////////////////////////////////
/// \brief exp
/// \param mu
/// \return
///
Eigen::Matrix4f exp(const Eigen::VectorXf& mu)
{

    static const float one_6th = 1.0/6.0;
    static const float one_20th = 1.0/20.0;

    Eigen::Vector3f w = mu.segment(3,3);
    float theta_sq = w.transpose()*w;
    float theta = sqrt(theta_sq);
    float A, B;


    Eigen::Matrix4f result = Eigen::Matrix4f::Identity();
    Eigen::Vector3f w2 = mu.segment(0,3);
    Eigen::Vector3f cross_v = w.cross(w2);


    if (theta_sq < 1e-8) {

        A = 1.0 - one_6th * theta_sq;
        B = 0.5;
        result.block(0,3,3,1) = mu.segment(0,3) + 0.5 * cross_v;
        //result.get_translation() = mu.template slice<0,3>() + 0.5 * cross;
    } else {
        float C;
        if (theta_sq < 1e-6) {
            C = one_6th*(1.0 - one_20th * theta_sq);
            A = 1.0 - theta_sq * C;
            B = 0.5 - 0.25 * one_6th * theta_sq;
        } else {
            const float inv_theta = 1.0/theta;
            A = sin(theta) * inv_theta;
            B = (1 - cos(theta)) * (inv_theta * inv_theta);
            C = (1 - A) * (inv_theta * inv_theta);
        }

//        std::cout << "la 1.5c " << " A "<< A << " B "<< B << " C " << C << std::endl;
        result.block(0,3,3,1) = mu.segment(0,3) + B * cross_v + C * w.cross(cross_v);

    }

    Eigen::Matrix3f R = result.block(0,0,3,3);
    rodrigues_so3_exp(w, A, B, R);
    result.block(0,0,3,3) = R;

    return result;
}


/////////////////////////////////////////////
/// \brief rodrigues_so3_exp
/// \param w
/// \param A
/// \param B
/// \param R
///
void rodrigues_so3_exp(const Eigen::Vector3f w, const float A, const float B, Eigen::Matrix3f &R)
{
    {
        const double wx2 = w[0]*w[0];
        const double wy2 = w[1]*w[1];
        const double wz2 = w[2]*w[2];

        R(0,0) = 1.0 - B*(wy2 + wz2);
        R(1,1) = 1.0 - B*(wx2 + wz2);
        R(2,2) = 1.0 - B*(wx2 + wy2);
    }
    {
        const double a = A*w[2];
        const double b = B*(w[0]*w[1]);
        R(0,1) = b - a;
        R(1,0) = b + a;
    }
    {
        const double a = A*w[1];
        const double b = B*(w[0]*w[2]);
        R(0,2) = b + a;
        R(2,0) = b - a;
    }
    {
        const double a = A*w[0];
        const double b = B*(w[1]*w[2]);
        R(1,2) = b - a;
        R(2,1) = b + a;
    }
}
