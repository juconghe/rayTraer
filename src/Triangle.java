// Triangle class
// defines a Triangle shape

import javax.vecmath.*;

public class Triangle extends Shape {
	private Vector3f p0, p1, p2;	// three vertices make a triangle
	private Vector3f n0, n1, n2;	// normal at each vertex

	public Triangle() {
	}
	public Triangle(Vector3f _p0, Vector3f _p1, Vector3f _p2, Material mat) {
		p0 = new Vector3f(_p0);
		p1 = new Vector3f(_p1);
		p2 = new Vector3f(_p2);
		material = mat;
		Vector3f normal = new Vector3f();
		Vector3f v1 = new Vector3f();
		Vector3f v2 = new Vector3f();
		v1.sub(p1, p0);
		v2.sub(p2, p0);
		normal.cross(v1, v2);
		normal.normalize();				// compute default normal:
		n0 = new Vector3f(normal);		// the normal of the plane defined by the triangle
		n1 = new Vector3f(normal);
		n2 = new Vector3f(normal);
	}
	public Triangle(Vector3f _p0, Vector3f _p1, Vector3f _p2,
					Vector3f _n0, Vector3f _n1, Vector3f _n2,
					Material mat) {
		p0 = new Vector3f(_p0);
		p1 = new Vector3f(_p1);
		p2 = new Vector3f(_p2);
		material = mat;
		n0 = new Vector3f(_n0);		// the normal of the plane defined by the triangle
		n1 = new Vector3f(_n1);
		n2 = new Vector3f(_n2);
	}
	public HitRecord hit(Ray ray, float tmin, float tmax) {

		/* YOUR WORK HERE: complete the triangle's intersection routine
		 * Normal should be computed by a bilinear interpolation from n0, n1 and n2
		 * using the barycentric coordinates: alpha, beta, (1.0 - alpha - beta) */
		Vector3f d = ray.getDirection();

		Matrix3f original = new Matrix3f(p1.x - p0.x, p2.x - p0.x, -d.x,
				                       p1.y - p0.y, p2.y - p0.y, -d.y,
										p1.z - p0.z, p2.z - p0.z, -d.z);

		Vector3f answer = new Vector3f(ray.getOrigin().x - p0.x, ray.getOrigin().y - p0.y,
										ray.getOrigin().z - p0.z);

        Matrix3f mx = new Matrix3f(p1.x - p0.x, p2.x - p0.x, -d.x,
                                    p1.y - p0.y, p2.y - p0.y, -d.y,
                                    p1.z - p0.z, p2.z - p0.z, -d.z);
        mx.setColumn(0, answer);
        Matrix3f my = new Matrix3f(p1.x - p0.x, p2.x - p0.x, -d.x,
                                    p1.y - p0.y, p2.y - p0.y, -d.y,
                                    p1.z - p0.z, p2.z - p0.z, -d.z);
        my.setColumn(1, answer);
        Matrix3f mz = new Matrix3f(p1.x - p0.x, p2.x - p0.x, -d.x,
                                    p1.y - p0.y, p2.y - p0.y, -d.y,
                                    p1.z - p0.z, p2.z - p0.z, -d.z);
        mz.setColumn(2, answer);

        float dOriginal = original.determinant();
        float beta = mx.determinant() / dOriginal;
        float gamma = my.determinant() / dOriginal;
        float alpha = 1 - beta - gamma;
        float t = mz.determinant() / dOriginal;

        if (beta < 0 || beta > 1 || gamma < 0 || gamma > 1 || alpha < 0 || alpha > 1) {
            return null;
        }

        if (t < tmin || t > tmax) {
            return null;
        }


        Vector3f normal = new Vector3f(alpha * n0.x + beta * n1.x + gamma * n2.x,
                                        alpha * n0.y + beta * n1.y + gamma * n2.y,
                                        alpha * n0.z + beta * n1.z + gamma * n2.z);
        HitRecord rec = new HitRecord();
        rec.pos = ray.pointAt(t);
        rec.t = t;
        rec.material = material;
        rec.normal = normal;
        rec.normal.normalize();

        return rec;
    }
}
