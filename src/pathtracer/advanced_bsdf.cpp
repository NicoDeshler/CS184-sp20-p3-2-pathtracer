#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>


using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Implement MirrorBSDF

    *pdf = 1.;
    reflect(wo, wi);
    return reflectance / abs_cos_theta(*wi);

}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
  // TODO: proj3-2, part 3
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
    return  exp(-1. * pow(tan(acos(h.z)) / alpha, 2)) / (PI * pow(alpha, 2) * pow(h.z, 4));
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
  // TODO: proj3-2, part 3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Spectrum.
    double cosTheta_i = abs_cos_theta(wi);
    Spectrum A = (eta * eta + k * k);
    Spectrum B = 2. * eta * cosTheta_i;
    double C = cosTheta_i * cosTheta_i;

    Spectrum Rs = (A - B + C) / (A + B + C);
    Spectrum Rp = (A * C - B + 1.) / (A * C + B + 1.);
    Spectrum F = (Rs + Rp) / 2.;

  return F;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO: proj3-2, part 3
  // Implement microfacet model here.
    if (wo.z > 0 && wi.z > 0) {
        Vector3D h = (wo + wi) / (wo + wi).norm();// half vector
        return F(wi) * G(wo, wi) * D(h) / 4. / wo.z / wi.z;
    }
    else {
        return Spectrum(0.);
    }
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: proj3-2, part 3
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.
  
  // Sampling half angle vector in spherial coordinates
    Vector2D r = sampler.get_sample();
    double r1 = r.x;
    double r2 = r.y;

    double theta_h = atan(sqrt(-1. * pow(alpha, 2) * log(1. - r1)));
    double phi_h = 2. * PI * r2;
    Vector3D h(sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h));
    
  // Assign wi based on the sampled half vector
    *wi = (2. * dot(h, wo)) * h - wo;
    wi->normalize();

  // Check sampled wi validity
    if ((*wi).z < 0) {
        *pdf = 0.;
        return Spectrum(0.);

    }
    
  // Determine pdf of having sampled the final wi
    double p_theta = (2. * sin(theta_h) / (pow(alpha, 2) * pow(cos(theta_h), 3))) * exp(-1. * pow(tan(theta_h)/alpha, 2));
    double p_phi = 1. / 2. / PI;
    double p_h = p_theta * p_phi / sin(theta_h);
    *pdf = p_h / 4. / dot(*wi, h);

  return MicrofacetBSDF::f(wo, *wi);
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Implement RefractionBSDF
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305
    bool TIR = !refract(wo, wi, ior);
    
    double wo_z = wo.z;
    double eta = (wo_z > 0) ? 1. / ior : ior;

    Spectrum out;

    if (TIR) {
        reflect(wo, wi);
        *pdf = 1;
        return reflectance / abs_cos_theta(*wi);
    }
    else {
        double R0 = pow((ior - 1) / (ior + 1), 2);
        double R = R0 + (1 - R0) * pow((1 - abs(wo_z)), 5); // Schlick's approximation for the Fresnel Factor
        if (coin_flip(R)) {
            reflect(wo, wi);
            *pdf = R;
            return R * reflectance / abs_cos_theta(*wi);
        }
        else {
            *pdf = 1 - R;
            return (1 - R) * transmittance / abs_cos_theta(*wi) * pow(eta,2);
        }
    }

  
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
  // TODO:
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
    Matrix3x3 ref(-1., 0., 0.,
        0., -1., 0.,
        0., 0., 1.);
    *wi = ref * wo;
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {
  // TODO:
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
    
    double wo_z = wo.z;

    // Define ratio of refractive indices based on entering/exiting refractive material
    float eta = (wo_z > 0) ? 1. / ior : ior;
    float wi_z2 = 1 - pow(eta, 2) * (1 - pow(wo_z, 2));
    bool TIR = wi_z2 < 0;   // total internal reflection condition

    if (!TIR) {
        wi->x = -1. * eta * wo.x;
        wi->y = -1. * eta * wo.y;
        wi->z = -1. * ((double)(wo_z > 0) - (double)(wo_z < 0)) * sqrt(wi_z2);
        wi->normalize();
        return true;
    } 
    return false;
 }

} // namespace CGL
