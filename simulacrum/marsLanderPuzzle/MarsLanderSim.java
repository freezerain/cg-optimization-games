package simulacrum.marsLanderPuzzle;

import src.DrawerController;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Scanner;

import static java.lang.Thread.sleep;

public class MarsLanderSim {
    private final static int MARS_WIDTH = 7000;
    private final static int MARS_HEIGHT = 3000;
    private final static int VERTICAL_SPEED_LIMIT = 40; // m/s
    private final static int HORIZONTAL_SPEED_LIMIT = 20; // m/s
    private final static File TEST_FILE = new File(MarsLanderSim.class.getResource("test1.txt").getPath());
    private static MarsLander ms;
    private static DrawerController dc;
    private final static double[] testLine = { 507.0, 653.0, 71.0, 337.0, 209.0, 119.0 };
    private final static double[] testLine2 = { 473.0, 671.0, 226.0, 131.0, 79.0 , 75.0};
    private final static double[] testLine3 = { 859.0, 699.0, 615.0, 505.0, 247.0, 136.0 };
    

    public static void main(String[] args) throws FileNotFoundException, InterruptedException {
        Scanner sc = new Scanner(TEST_FILE);
        double[] landscape = Arrays.stream(sc.nextLine().split("(, )"))
                .mapToDouble(Double::parseDouble)
                .toArray();
        int[] s1 = Arrays.stream(sc.nextLine().split(" ")).mapToInt(Integer::parseInt).toArray();
        ms = new MarsLander(s1[0],s1[1],s1[2],s1[3],s1[4],s1[5],s1[6]);
        dc = new DrawerController();
        dc.startApp(landscape);

        for (int i = 0; i < 50; i++){
           nextTurn();
           sleep(500);
        }
    }
    public static void nextTurn(){
        ms.evaluateTurn(-20,3);
        dc.updateCoord(ms.x, ms.y, ms.hSpeed, ms.vSpeed, ms.rotateAngle, ms.power, ms.fuel, new double[][]{testLine, testLine2, testLine3});
    }

}

