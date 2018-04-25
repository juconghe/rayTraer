// Sphere class
// defines a Sphere shape

import javax.vecmath.*;

public class Sphere extends Shape {
	private Vector3f center;	// center of sphere
	private float radius;		// radius of sphere

	public Sphere() {
	}
	public Sphere(Vector3f pos, float r, Material mat) {
		center = new Vector3f(pos);
		radius = r;
		material = mat;
	}
	public HitRecord hit(Ray ray, float tmin, float tmax) {

		/* YOUR WORK HERE: complete the sphere's intersection routine */
		Vector3f oriToC = new Vector3f(ray.getOrigin());
        Vector3f direction = new Vector3f(ray.getDirection());

		oriToC.sub(center);

		float a = direction.lengthSquared();
		float b = 2* oriToC.dot(direction);
		float c = oriToC.lengthSquared() - (radius * radius);

		float rootValue = (float) (b * b - 4.0 * a * c);
		if (rootValue <= 0) {
		    return null;
        }

		float sqRootValue = (float) Math.sqrt(rootValue);
		float t1 = (-b + sqRootValue) / (2.0f * a);
		float t2 = (-b - sqRootValue) / (2.0f * a);

		float t;
        if (tmin<=t1 && t1<=tmax) {
            if (tmin<=t2 && t2<=tmax && t2<t1){
                t=t2;
            }else {
                t = t1;
            }
        } else if (tmin<=t2 && t2<=tmax) {
            t = t2;
        } else {
            return null;
        }

		Vector3f normal = new Vector3f(ray.pointAt(t));
		normal.sub(center);
		HitRecord rec = new HitRecord();
		rec.pos = ray.pointAt(t);
		rec.t = t;
		rec.material = material;
		rec.normal = normal;
		rec.normal.normalize();

		return rec;
	}
}
