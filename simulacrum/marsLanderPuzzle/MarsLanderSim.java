package simulacrum.marsLanderPuzzle;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Scanner;

class MarsLanderSim {
    private final static int MARS_WIDTH = 7000;
    private final static int MARS_HEIGHT = 3000;
    private final static int VERTICAL_SPEED_LIMIT = 40; // m/s
    private final static int HORIZONTAL_SPEED_LIMIT = 20; // m/s
    private final static File TEST_FILE = new File(MarsLanderSim.class.getResource("test1.txt").getPath());

    public static void main(String[] args) throws FileNotFoundException {
        Scanner sc = new Scanner(TEST_FILE);
        int[] landscape = Arrays.stream(sc.nextLine().split("[ \n]"))
                .mapToInt(Integer::parseInt)
                .toArray();
        int[] s1 = Arrays.stream(sc.nextLine().split(" ")).mapToInt(Integer::parseInt).toArray();
        MarsLander ms = new MarsLander(new GameState(s1[0],s1[1],s1[2],s1[3],s1[4],s1[5],s1[6]));
        for (int i = 0; i < 6; i++){
            ms.evaluateTurn(-20,3);
        }
    }

}

