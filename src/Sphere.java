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
		Vector3f oriToC = ray.getOrigin();
		oriToC.sub(center);

		float a = ray.getDirection().lengthSquared();
		float b = 2 * oriToC.dot(ray.getDirection());
		float c = oriToC.lengthSquared() - (radius * radius);

		float sqRootValue = (float) Math.sqrt(Math.pow(b, 2) - 4.0 * a * c);
		float t1 = (-b + sqRootValue) / (2.0f * a);
		float t2 = (-b - sqRootValue) / (2.0f * a);

		float t;
		if (tmin <= t1 && t1 <= tmax) {
			t = t1;
		} else if (tmin <= t2 && t2 <= tmax) {
			t = t2;
		} else {
			return null;
		}

		Vector3f interPoint = ray.pointAt(t);
		interPoint.sub(center);
		HitRecord rec = new HitRecord();
		rec.pos = ray.pointAt(t);
		rec.t = t;
		rec.material = material;
		rec.normal = interPoint;
		rec.normal.normalize();

		return rec;
	}
}
