package com.expression.parser;

import com.expression.parser.util.Point;
import org.junit.runners.ParentRunner;

import java.util.Arrays;
import java.util.Scanner;

public class Main {

    static Scanner sc = new Scanner(System.in);
    public static void main(String[] args) {
        int n = 2;
        Point[] x_k = new Point[2];
        x_k[0] = new Point( "x", 0.0); // початкові значення
        x_k[1] = new Point("y", 1.9);

        String[] vars = new String[n];
        vars[0] = x_k[0].getVar();
        vars[1] = x_k[1].getVar();
        methodNewton(x_k, vars);
    }

    public static void print(double[] arr){
        for(double item : arr){
            System.out.print(item + " ");
        }
        System.out.println();
    }

    public static void print(double[][] matrix){
        for(int i = 0; i < matrix.length; i++){
            for(int j = 0; j < matrix.length; j++){
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();
    }

    public static boolean isEpsilon(Point[] x_k, double[] k_p, double e){
        for(int i = 0; i < x_k.length; i++){
            if(x_k[i].getValue() - k_p[i] > e)
                return false;
        }
        return true;
    }

    // TODO: this method do some magic lol
    public static void methodNewton(Point[] x_k, String[] vars){
        int iter = 0;
        double e = 0.0001;
        double[] F = new double[x_k.length];
        double[] k_p = new double[x_k.length];
        String[] str_f = new String[x_k.length];
        str_f[0] = "sin(x)+sqrt(2*y^3)-4";
        str_f[1] = "tan(x)-y^2+4";

        do{
            for(int i = 0; i < F.length; i++){
                //str_f[i] = sc.next();
                F[i] = Parser.eval(str_f[i], x_k).getValue() * -1;
            }
            System.out.println("F:");
            print(F);

            double[][] J = Jacobian(x_k, str_f, vars);
            System.out.println("J:");
            print(J);

            System.out.println("GaussEl:");
            double[] delta_k = methodGaussEl(J, F);
            System.out.println("delta:");
            print(delta_k);

            System.out.println("points:");

            for(int i = 0; i < k_p.length; i++){
                k_p[i] = x_k[i].getValue();
            }
            for(int i = 0; i < x_k.length; i++){
                x_k[i].setValue(x_k[i].getValue() + delta_k[i]);
                System.out.println(x_k[i].getValue());
            }
            System.out.println("-------------------------------------");
            iter++;

        }while(!isEpsilon(x_k, k_p, e));
        System.out.println("iterations: " + iter);

        System.out.println("Check results: ");
        double rez_f1 = Parser.eval(str_f[0], x_k).getValue();
        double rez_f2 = Parser.eval(str_f[1], x_k).getValue();

        System.out.println(rez_f1 + "\n" + rez_f2);
    }

    public static double[] methodGaussEl(double[][] matrix, double[] b){
        for(int k = 0; k < matrix.length - 1; k++){

            for(int i = k + 1; i < matrix.length; i++){

                double a = (matrix[i][k] / matrix[k][k]) * -1;
                b[i] += b[k] * a;
                for(int j = 0; j < matrix.length; j++){
                    matrix[i][j] += (matrix[k][j] * a);
                }
            }
        }
        print(matrix);
        System.out.println("b:");
        print(b);

        int n = matrix.length - 1;
        double[] delta = new double[b.length];
        Arrays.fill(delta,1);
        for(int i = n; i >= 0; i--){

            double Sum = 0;
            for(int j = n; j > i; j--){
                Sum += matrix[i][j] * delta[j];
            }
            delta[i] = (b[i] - Sum) / matrix[i][i];
        }
        return delta;
    }

    // TODO: this method return a matrix with derivative df of k
    public static double[][] Jacobian(Point[] x_k, String[] str_f, String[] vars){

        double[][] J = new double[x_k.length][x_k.length];

        Double[] val = new Double[x_k.length];
        Arrays.fill(val,1.0);

        for(int i = 0; i < J.length; i++){

            val[i] = x_k[i].getValue();
            for(int j = 0; j < J.length; j++){
                double f = Parser.eval(str_f[j], vars, val);
                val[i] += 0.05;
                J[j][i] = (float)((Parser.eval(str_f[j], vars, val) - f) / 0.05);
            }
            val[i] = 1.0;
        }
        return J;
    }
}

