import java.util.*;
import java.util.stream.Collectors;

class Player {
    static final int[] FIBONACCI = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
            1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811,
            514229, 832040, 1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817,
            39088169, 63245986, 102334155};

    public static void main(String args[]) {
        Scanner in = new Scanner(System.in);
        List<Chromosome> population = new ArrayList<>();
        while (true) {
            Timer timer = new Timer();
            Actor player = new Actor(-1, in.nextInt(), in.nextInt());
            List<Actor> humans = new ArrayList<>();
            List<Actor> zombies = new ArrayList<>();
            int humanCount = in.nextInt();
            for (int i = 0; i < humanCount; i++)
                humans.add(new Actor(in.nextInt(), in.nextInt(), in.nextInt()));
            int zombieCount = in.nextInt();
            for (int i = 0; i < zombieCount; i++){
                zombies.add(new Actor(in.nextInt(), in.nextInt(), in.nextInt()));
                in.nextInt();
                in.nextInt();
            }
            System.err.println("before eval:" + timer);
            List<Chromosome> generated = findMove(timer, player, humans, zombies);
            generated.sort(Comparator.reverseOrder());
            Genome best = generated.get(0).gene.get(0);
            System.out.println(best.x +" "+ best.y);
            System.err.println("after eval:" + timer);
        }

    }
    private static List<Chromosome> mutate(Timer t,  Actor player, List<Actor> humans, List<Actor> zombies, List<Chromosome> population){
        int counter = 0;
        List<Chromosome> result = new ArrayList<>();
        for (int i = 0; i < population.size() && t.getTime()<90; i++){
            
        }
        return result;
    }
    private static List<Chromosome> findMove(Timer t, Actor player, List<Actor> humans, List<Actor> zombies) {
        int counter = 0;
        //Generate random
        List<Chromosome> population = new ArrayList<>();
        while (t.getTime() < 90 && counter<100000) {
            Actor copyPlayer = new Actor(player);
            List<Actor> copyHumans = humans.stream()
                    .map(Actor::new)
                    .collect(Collectors.toList());
            List<Actor> copyZombies = humans.stream()
                    .map(Actor::new)
                    .collect(Collectors.toList());
            Evaluator evaluator = new Evaluator(copyPlayer, copyHumans, copyZombies, 100);
            Chromosome eval = evaluator.eval();
            if(eval.humansAlive>0) {
                population.add(eval);
                System.out.println("find 1 solution");
            }
            counter++;
        }
        System.err.println("evaluated:"+counter + " before:" + t.getTime());
        return population;
    }

    static class Evaluator {
        Actor player;
        List<Actor> humans;
        List<Actor> zombies;
        int maxTurn;
        List<Genome> genes = new ArrayList<>();
        Random r = new Random();
        int score = 0;

        public Evaluator(Actor player, List<Actor> humans, List<Actor> zombies, int maxTurn) {
            this.player = player;
            this.humans = humans;
            this.zombies = zombies;
            this.maxTurn = maxTurn;
        }

        public Chromosome eval() {
            for (int i = 0; i < maxTurn; i++){
                madeTurn();
                if (zombies.isEmpty() || humans.isEmpty()) {
                    break;
                }
            }
            return new Chromosome(genes, score, humans.size(),0);
        }
        
        private void madeTurn() {
            int playerNextX = r.nextInt(16000);
            int playerNextY = r.nextInt(9000);
            genes.add(new Genome(playerNextX, playerNextY));
            //1) Move Zombies
            moveZombies();
            //2) Move Player
            player.moveToNext(playerNextX, playerNextY, 1000);
            //3) Player kills zombies
            killZombies();
            //4) Zombies eat humans
            killHumans();
        }

        private void killHumans() {
            for (Actor z: zombies){
                z.getInRange(humans, 0).forEach(h -> humans.remove(h));
            }
        }

        private void killZombies() {
            List<Actor> inRange = player.getInRange(zombies, 2000);
            if (inRange.size() == 0) return;
            int fibSum = 0;
            for (int i = 0; i < inRange.size(); i++){
                fibSum += FIBONACCI[i];
            }
            score += 10 * (humans.size() * humans.size()) * fibSum;
            inRange.forEach(z -> zombies.remove(z));
        }

        private void moveZombies() {
            for (Actor z: zombies){
                Actor closestHuman = z.getClosest(humans);
                Actor closest = closestHuman != null &&
                                z.getDistanceSqrt(closestHuman) < z.getDistanceSqrt(player) ?
                        closestHuman : player;
                z.moveToNext(closest.x, closest.y, 400);
            }
        }
    }

    static class Chromosome implements Comparable<Chromosome>{
        List<Genome> gene;
        int score;
        int humansAlive;
        int canBeSaved;

        public Chromosome(List<Genome> gene, int score, int humansAlive, int canBeSaved) {
            this.gene = gene;
            this.score = score;
            this.humansAlive = humansAlive;
            this.canBeSaved = canBeSaved;
        }

        public Chromosome compare(Chromosome candidate) {
            if(score != candidate.score) return score>candidate.score? this:candidate;
            if(humansAlive!=candidate.humansAlive) return humansAlive>candidate.humansAlive? this:candidate;
            if(canBeSaved!=candidate.canBeSaved) return canBeSaved>candidate.canBeSaved? this:candidate;
            if(gene.size()!=candidate.gene.size()) return gene.size()<candidate.gene.size()? this:candidate;
            return this;
        }

        @Override
        public int compareTo(Chromosome o) {
            if(humansAlive!=o.humansAlive) return Integer.compare(humansAlive, o.humansAlive);
            if(score != o.score) return Integer.compare(score, o.score);
            if(gene.size()!=o.gene.size()) return Integer.compare(o.gene.size(),gene.size());
            return 0;
        }
    }

    static class Genome {
        int x;
        int y;

        public Genome(int x, int y) {
            this.x = x;
            this.y = y;
        }
    }

    static class Actor {
        int id;
        int x;
        int y;

        public Actor(int id, int x, int y) {
            this.id = id;
            this.x = x;
            this.y = y;
        }
        public Actor(Actor a){
            this.id = a.id;
            this.x=a.x;
            this.y=a.y;
        }

        public int getDistanceSqrt(Actor target) {
            return (target.x - x) * (target.x - x) + (target.y - y) * (target.y - y);
        }

        public int getDistanceSqrt(int targetX, int targetY) {
            return (targetX - x) * (targetX - x) + (targetY - y) * (targetY - y);
        }

        public List<Actor> getInRange(List<Actor> actorList, int range) {
            return actorList.stream()
                    .filter(a -> this.getDistanceSqrt(a) <= (range * range))
                    .collect(Collectors.toList());
        }

        public Actor getClosest(List<Actor> actorList) {
            return actorList.stream()
                    .min(Comparator.comparingInt(this::getDistanceSqrt))
                    .orElse(null);
        }


        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Actor actor = (Actor) o;
            return id == actor.id;
        }

        @Override
        public int hashCode() {
            return id;
        }

        public void moveToNext(int nextX, int nextY, int movementDistance) {
            if (getDistanceSqrt(nextX, nextY) <= (movementDistance * movementDistance)) {
                x = nextX;
                y = nextY;
            } else {
                double vectorX = nextX - x;
                double vectorY = nextY - y;
                double normaliser = Math.max(Math.abs(vectorX), Math.abs(vectorY));
                x = (int) (x + (vectorX / normaliser * movementDistance));
                y = (int) (y + (vectorY / normaliser * movementDistance));
            }
        }
    }

    static class Timer {
        long start;

        public Timer() {
            start = System.currentTimeMillis();
        }

        public long getTime() {
            return System.currentTimeMillis() - start;
        }

        @Override
        public String toString() {
            return (System.currentTimeMillis() - start) + "ms.";
        }
    }
}