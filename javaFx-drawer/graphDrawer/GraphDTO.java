package graphDrawer;

public class GraphDTO {
    String name;
    boolean isDouble;
    boolean boolVar;
    double doubleVar;
    
    int time;
    int firstSum;
    int tenSum;
    double bestFit;
    double lastFit;
    int failed;

    public GraphDTO(String name, boolean isDouble, boolean boolVar, double doubleVar, int time,
                    int firstSum, int tenSum, double bestFit, double lastFit, int failed) {
        this.name = name;
        this.isDouble = isDouble;
        this.boolVar = boolVar;
        this.doubleVar = doubleVar;
        this.time = time;
        this.firstSum = firstSum;
        this.tenSum = tenSum;
        this.bestFit = bestFit;
        this.lastFit = lastFit;
        this.failed = failed;
    }
}
