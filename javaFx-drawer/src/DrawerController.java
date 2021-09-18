package src;

import graphDrawer.GraphDTO;
import graphDrawer.GraphDrawerApp;
import javafx.application.Application;
import marsLander.LanderPath;
import marsLander.MarsLanderApp;
import java.util.List;
import java.util.Map;

public class DrawerController {
    private MarsLanderApp app;

    public void updateCoord(int gen, List<LanderPath> paths) {
        app.updateObservable(gen, paths);
    }
    
    public void showGraph(String testName, Map<String, Map<Boolean, Number>> data){
        startGraphApp().showBarChart(testName, data);
    }

    public void showLineChart(String testName, Map<String, Map<Number, Number>> data){
        startGraphApp().showLineChart(testName, data);
    }
    public void showMultiLineChart(String testName, Map<String, Map<double[], Number>> data){
        startGraphApp().showMultiLineChart(testName, data);
    }
    
    private GraphDrawerApp startGraphApp(){
        new Thread(() -> Application.launch(GraphDrawerApp.class)).start();
        return GraphDrawerApp.waitForStart();
    }

    public void startApp(double[] landscape) {
        long time = System.currentTimeMillis();
        new Thread(() -> Application.launch(MarsLanderApp.class)).start();
        double[] normLand = new double[landscape.length];
        for (int i = 0; i < landscape.length; i++){
            if (i % 2 != 0) normLand[i] = 750.0 - landscape[i] / 4.0;
            else normLand[i] = landscape[i] / 4.0;
        }
        app = MarsLanderApp.waitForStart();
        app.buildLandscape(normLand);
        System.out.println("await finished: " + (System.currentTimeMillis() - time) + "ms.");
    }

}
