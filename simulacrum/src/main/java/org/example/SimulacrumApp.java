package org.example;

import javafx.application.Application;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Label;
import javafx.scene.control.Separator;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.paint.Paint;
import javafx.scene.shape.*;
import javafx.stage.Stage;
import javafx.stage.StageStyle;

import java.util.List;

public class SimulacrumApp extends Application {

    @Override
    public void start(Stage stage) throws Exception {
        stage.setTitle("Simulacrum");
        stage.initStyle(StageStyle.UNIFIED);
        
        Pane pane = new Pane();
        BackgroundFill b1 = new BackgroundFill(Color.color(0.0, 0.0, 0.0), new CornerRadii(0),
                new Insets(0));
        BackgroundFill b2 = new BackgroundFill(Color.color(0.5, 0.5, 0.5), new CornerRadii(5),
                new Insets(0));
        pane.setBackground(new Background(b1,b2));
        
        Circle circle = new Circle();
        circle.setRadius(20);
        circle.setLayoutX(700);
        circle.setLayoutY(200);
        //pane.getChildren().add(vbox);
        pane.getChildren().add(circle);
        Canvas canvas = new Canvas();
        canvas.prefWidth(Double.MAX_VALUE);
        canvas.prefHeight(Double.MAX_VALUE);
        canvas.setWidth(1280);
        canvas.setHeight(720);
        GraphicsContext gc = canvas.getGraphicsContext2D();
        gc.setFill(Color.color(0,0,0));
        
        gc.setLineWidth(1.0);
        for (double i = 0.5; i < 200; i+=10){
            gc.moveTo(i, 0);
            gc.lineTo(i, 600);
            gc.stroke();
        }
        for (double i = 0.5; i < 200; i+=10){
            gc.moveTo(0, i);
            gc.lineTo(600, i);
            gc.stroke();
        }
        Line line = new Line(500,200,750,600);
        line.setStrokeWidth(5);
        pane.getChildren().add(line);
        Line line2 = new Line(450,150,700,550);
        line2.setStrokeWidth(10);
        pane.getChildren().add(line2);
        //gc.strokeLine(500, 200 , 750 , 600);
        canvas.setHeight(220);
        canvas.setWidth(1800);
        pane.getChildren().add(canvas);
        stage.setWidth(1800);
        stage.setHeight(1000);
        stage.setScene(new Scene(pane));
        stage.show();
    }

    public static void main(String[] args) {
        launch();
    }
}
