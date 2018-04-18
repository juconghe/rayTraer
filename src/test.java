import javax.vecmath.*;

public class test {
    public static void main(String[] args) {
        Matrix3f original = new Matrix3f(1,2,1,3,2,0,0,1,0);
        Vector3f answer = new Vector3f(1,2,3);

        Matrix3f mx = new Matrix3f(1,2,1,3,2,0,0,1,0);
        mx.setColumn(0, answer);
        Matrix3f my = new Matrix3f(1,2,1,3,2,0,0,1,0);
        my.setColumn(1, answer);
        Matrix3f mz = new Matrix3f(1,2,1,3,2,0,0,1,0);
        mz.setColumn(2, answer);

        float dOriginal = original.determinant();
        float dx = mx.determinant();
        float dy = my.determinant();
        float dz = mz.determinant();

        System.out.println("x: " + dx/dOriginal);
        System.out.println("Y: " + dy/dOriginal);
        System.out.println("Z: " + dz/dOriginal);

    }
}
