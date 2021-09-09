package src;

import javafx.application.Application;
import marsLander.MarsLanderApp;

import java.lang.reflect.Array;

public class DrawerController {
    int maxWidth = 1800;
    int maxHeight = 1000;
    private MarsLanderApp app;
    int normalizedWidth = 1750;
    int normalizedHeight = 750;
    
    public void updateCoord(double x, double y, double hSpeed, double vSpeed, int angle, int power, int fuel, double[][] paths){
        int normX = (int) Math.round(x / 4.0);
        int normY = normalizedHeight - (int) Math.round(y / 4.0);
        app.updateObservable(normX, normY,  (int) Math.round(hSpeed),  (int) Math.round(vSpeed), angle, power, fuel, paths);
    }
    
    public void startApp(double[] landscape){
        long time = System.currentTimeMillis();
        new Thread(() -> Application.launch(MarsLanderApp.class)).start();
        app = MarsLanderApp.waitForStart();
        app.buildLandscape(landscape);
        System.out.println("await finished: " + (System.currentTimeMillis() - time) + "ms.");
    }
    
}
