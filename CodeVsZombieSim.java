import marsLander.Lander;
import marsLander.LanderPath;
import src.DrawerController;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static java.lang.Thread.sleep;

public class CodeVsZombieSim {
    static int[] player = {500, 4500};
    static int[] humans = {100,4000,130,5000,10,4500,500,3500,10,5500,100,3000};
    static int[] zombies = {8000,4500,9000,4500,10000,4500,11000,4500,12000,4500,
            13000,4500,14000,4500,15000,3500,14500,2500,15900,500};

    public static void main(String[] args) throws InterruptedException {
        evaluate(player,humans, zombies);
    }

    private static void evaluate(int[] player, int[] humans, int[] zombies) throws InterruptedException {
        Player.Genetic.Settings settings = new Player.Genetic.Settings();
        Player.Genetic.Solver solver = new Player.Genetic.Solver(settings);
        List<Player.Genetic.Individual> pop = new ArrayList<>();
        int id = 0;
        List<Player.Actor> humansActor = new ArrayList<>();
        for (int i = 0; i < humans.length; i+=2){
            humansActor.add(new Player.Actor(id++, humans[i], humans[i+1]));
        }
        List<Player.Actor> zombiesActor = new ArrayList<>();
        for (int i = 0; i < zombies.length; i+=2){
            zombiesActor.add(new Player.Actor(id++, zombies[i], zombies[i+1]));
        }
        Player.Actor playerActor = new Player.Actor(-1, player[0], player[1]);
        DrawerController dc = new DrawerController();
        Player.Genetic.State startState = new Player.Genetic.State(playerActor, humansActor, zombiesActor);
        while((!startState.zombies.isEmpty()) && (!startState.humans.isEmpty())){
            int[] counter = new int[1];
         
            pop = solver.evolve(pop, startState, System.currentTimeMillis(),counter);
            System.err.println("crossovers: "+ counter[0]);
            var scores = new double[pop.size()];
            var paths = new int[pop.size()][];
            
            
            for (int i = 0; i < pop.size(); i++){
                Player.Genetic.Individual individual = pop.get(i);
                scores[i] = individual.state.humans.isEmpty() ? -1.0:     individual.fitness;
                List<Player.Genetic.Gene> genes = individual.genes;
                int[] path = new int[genes.size()*2];
                for (int j = 0; j < genes.size(); j++){
                    Player.Genetic.Gene gene = genes.get(j);
                    Player.Genetic.State state = individual.states.get(j);
                    int nextX = state.player.x;
                    int nextY = state.player.y;
                    path[j*2] = nextX;
                    path[j*2+1] = nextY;
                }
                paths[i] = path;
            }
            
            
            
            startState = pop.get(0).states.get(0);
            System.err.println("score: " + startState.score);
            humans = new int[startState.humans.size()*2];
            for (int i = 0; i < startState.humans.size(); i++){
                humans[i*2] = startState.humans.get(i).x;
                humans[i*2+1] = startState.humans.get(i).y;
            }
            zombies = new int[startState.zombies.size()*2];
            for (int i = 0; i < startState.zombies.size(); i++){
                zombies[i*2] = startState.zombies.get(i).x;
                zombies[i*2+1] = startState.zombies.get(i).y;
            }
            
            dc.updateCodeVsZombieApp(new int[]{startState.player.x, startState.player.y}, 
                    humans, zombies,paths,scores);
            sleep(1000);
        }
    }
    
}
