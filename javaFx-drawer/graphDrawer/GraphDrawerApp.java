package graphDrawer;

import javafx.application.Application;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.Text;
import javafx.stage.Stage;
import marsLander.CanvasGridGenerator;
import marsLander.MarsLanderApp;

import java.util.concurrent.CountDownLatch;

public class GraphDrawerApp extends Application {
    private static final CountDownLatch latch = new CountDownLatch(1);
    public static GraphDrawerApp appReference = null;
    private Pane pane = new Pane();
    
    public GraphDrawerApp() {
        appReference = this;
    }

    public static GraphDrawerApp waitForStart() {
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return appReference;
    }
    
    
    @Override
    public void start(Stage stage) throws Exception {
        stage.setTitle("Graph drawer");
        pane = new Pane();
        pane.setBackground(new Background(
                new BackgroundFill(Color.color(0.0, 0.0, 0.0), new CornerRadii(0), new Insets(0))));
        
        
        
        stage.setWidth(1280);
        stage.setHeight(720);
        stage.setScene(new Scene(pane));
        stage.show();
        latch.countDown();
    }
}
