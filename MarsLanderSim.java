import simulacrum.marsLanderPuzzle.MarsLander;
import src.DrawerController;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import static java.lang.Thread.sleep;

public class MarsLanderSim {
    private final static int MARS_WIDTH = 7000;
    private final static int MARS_HEIGHT = 3000;
    private final static int VERTICAL_SPEED_LIMIT = 40; // m/s
    private final static int HORIZONTAL_SPEED_LIMIT = 20; // m/s
    private final static File TEST_FILE = new File(MarsLanderSim.class.getResource("simulacrum" +
                                                                                   "/marsLanderPuzzle/test1.txt").getPath());
    private static MarsLander ms;
    private static DrawerController dc;
    private static final int TEST_NUMBER = 1;
    

    public static void main(String[] args) throws FileNotFoundException, InterruptedException {
        Scanner sc = new Scanner(TEST_FILE);
        for (int i = 0; i < TEST_NUMBER*2; i++){
            sc.nextLine();
        }
        double[] landscape = Arrays.stream(sc.nextLine().split(" "))
                .mapToDouble(Double::parseDouble)
                .toArray();
        int[] s1 = Arrays.stream(sc.nextLine().split(" ")).mapToInt(Integer::parseInt).toArray();
        Player.GameState gs = new Player.GameState(s1[0],s1[1],s1[2],s1[3],s1[4],s1[5],s1[6]);
        dc = new DrawerController();
        Player player = new Player();
        player.simulationInit(landscape);
        dc.startApp(landscape);
        sleep(3000);
        evaluateStatic(gs,player);
        //evaluateAndSimulate(gs, player);
        System.out.println("Algho completed!");
    }

    private static void evaluateAndSimulate(Player.GameState gs, Player player) throws InterruptedException {
        int counter = 0;
        while(!gs.isLanded){
            List<Player.Individual> indList = player.simulationStep(gs);
            if(indList.isEmpty()) {
                System.out.println("population is empty");
                break;
            }
            gs = indList.get(0).gameStateList.get(0);
            HashMap<double[], Double> map = new HashMap<>();
            for (int i = 0; i < indList.size(); i++){
                Player.Individual ind = indList.get(i);
                List<Player.GameState> gameState = ind.gameStateList;
                List<Double> newPath = new ArrayList<>();
                for (int i1 = 0; i1 < gameState.size(); i1++){
                    Player.GameState step = gameState.get(i1);
                    newPath.add(step.x);
                    newPath.add(step.y);
                }
                double[] pathCoord = newPath.stream().mapToDouble(d -> d).toArray();
                map.put(pathCoord, ind.fitnessScore);
            }
            
            dc.updateCoord(counter++, gs.x, gs.y, gs.hSpeed, gs.vSpeed, gs.angle, gs.power, gs.fuel, map);
            sleep(1000);
        }
    }

    private static void evaluateStatic(Player.GameState gs, Player player) throws InterruptedException {
        int counter = 0;
        while(true){
            List<Player.Individual> indList = player.simulationStep(gs);
            if(indList.isEmpty()) {
                System.out.println("population is empty");
                break;
            }
            HashMap<double[], Double> map = new HashMap<>();
            for (int i = 0; i < indList.size(); i++){
                Player.Individual ind = indList.get(i);
                List<Player.GameState> gameState = ind.gameStateList;
                List<Double> newPath = new ArrayList<>();
                for (int i1 = 0; i1 < gameState.size(); i1++){
                    Player.GameState step = gameState.get(i1);
                    newPath.add(step.x);
                    newPath.add(step.y);
                }
                double[] pathCoord = newPath.stream().mapToDouble(d -> d).toArray();
                map.put(pathCoord, ind.fitnessScore);
            }

            dc.updateCoord(counter++, gs.x, gs.y, gs.hSpeed, gs.vSpeed, gs.angle, gs.power, gs.fuel, map);
            sleep(300);
        }
    }
}

