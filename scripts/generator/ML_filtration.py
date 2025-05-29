import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from PROBESt.AI import LogisticRegressionModel, PerceptronModel, DeepNeuralNetworkModel
from PROBESt.filtration import train_filtration_AI, validate_filtration_AI, apply_filtration_AI

def main():
    # Load data
    data_path = 'data/databases/open/test_ML_database.csv'
    data = pd.read_csv(data_path)
    
    # Convert boolean 'type' column to numeric
    data['type'] = data['type'].astype(int)
    
    # Drop character columns but keep the target column
    numeric_columns = data.select_dtypes(include=[np.number]).columns
    data = data[numeric_columns]
    
    # Drop rows with NaN values
    initial_size = len(data)
    data = data.dropna()
    print(f"Dropped {initial_size - len(data)} rows with NaN values")
    
    # Split data into train, validation, and test sets
    train_data, temp_data = train_test_split(data, test_size=0.3, random_state=42)
    val_data, test_data = train_test_split(temp_data, test_size=0.5, random_state=42)
    
    print(f"Training set size: {len(train_data)}")
    print(f"Validation set size: {len(val_data)}")
    print(f"Test set size: {len(test_data)}")
    
    # Train and evaluate deep neural network model
    print("\nTraining Deep Neural Network model...")
    dnn_model = DeepNeuralNetworkModel(input_size=train_data.shape[1] - 1)
    trained_dnn_model, dnn_metrics = train_filtration_AI(dnn_model, train_data)
    
    print("\nDeep Neural Network metrics:")
    for metric, value in dnn_metrics.items():
        print(f"{metric}: {value:.4f}")
    
    # Validate deep neural network model
    print("\nValidating Deep Neural Network model...")
    dnn_val_metrics = validate_filtration_AI(trained_dnn_model, val_data, output_name='DNN.png')
    
    print("\nDeep Neural Network validation metrics:")
    for metric, value in dnn_val_metrics.items():
        print(f"{metric}: {value:.4f}")
    
    # Train and evaluate logistic regression model
    print("\nTraining Logistic Regression model...")
    lr_model = LogisticRegressionModel()
    trained_lr_model, lr_metrics = train_filtration_AI(lr_model, train_data)
    
    print("\nLogistic Regression metrics:")
    for metric, value in lr_metrics.items():
        print(f"{metric}: {value:.4f}")
    
    # Validate logistic regression model
    print("\nValidating Logistic Regression model...")
    lr_val_metrics = validate_filtration_AI(trained_lr_model, val_data, output_name='LR.png')
    
    print("\nLogistic Regression validation metrics:")
    for metric, value in lr_val_metrics.items():
        print(f"{metric}: {value:.4f}")
    
    # Apply best model to test set
    print("\nApplying best model to test set...")
    models = {
        "Deep Neural Network": (trained_dnn_model, dnn_val_metrics),
        "Logistic Regression": (trained_lr_model, lr_val_metrics)
    }
    
    best_model_name = max(models.keys(), key=lambda k: models[k][1]['f1'])
    best_model = models[best_model_name][0]
    
    print(f"\nUsing {best_model_name} model for final predictions")    
    test_predictions = apply_filtration_AI(best_model, test_data)
    
    # Save predictions
    output_dir = 'tests_outs'
    os.makedirs(output_dir, exist_ok=True)
    test_predictions.to_csv(os.path.join(output_dir, 'test_predictions.csv'), index=False)
    print(f"\nPredictions saved to {os.path.join(output_dir, 'test_predictions.csv')}")

if __name__ == '__main__':
    main() 