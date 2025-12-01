import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from PROBESt.AI import LogisticRegressionModel, DeepNeuralNetworkModel, TorchClassifier
from PROBESt.filtration import (
    train_filtration_AI, validate_filtration_AI, apply_filtration_AI,
    plot_combined_roc_curves, plot_combined_metrics, plot_learning_curves
)
from PROBESt.models_registry import (
    ShallowNet, WideNet, ResidualNet, GAILDiscriminator, TabTransformer,
    GAILDeep, GAILWide, GAILNarrow, GAILWithDropout
)

MODELS = {
    "ShallowNet": lambda n: TorchClassifier(ShallowNet(n), weight_pos=5),
    "WideNet": lambda n: TorchClassifier(WideNet(n), weight_pos=5),
    "ResidualNet": lambda n: TorchClassifier(ResidualNet(n), weight_pos=5),
    "GAIL": lambda n: TorchClassifier(GAILDiscriminator(n), weight_pos=5),
    "TabTransformer": lambda n: TorchClassifier(TabTransformer(n), weight_pos=5),
}

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
    
    # Get input size
    input_size = train_data.shape[1] - 1
    
    # Train and validate all models from MODELS dictionary
    results = {}
    trained_models = {}
    models_for_plots = {}  # For combined plots: (model, val_data)
    
    for name, constructor in MODELS.items():
        print(f"\n===== Training {name} =====")
        model = constructor(input_size)
        
        # For GAIL, track learning curves during training
        if name == "GAIL":
            X_train = train_data.drop(columns=['type'])
            y_train = train_data['type']
            model.train(X_train, y_train, epochs=100, batch_size=32, 
                       val_data=val_data, track_curves=True)
            trained = model
            # Get training metrics from validation (we'll use val_metrics for results)
            metrics = {}
        else:
            trained, metrics = train_filtration_AI(model, train_data)
        
        val_metrics = validate_filtration_AI(trained, val_data, output_name=f"{name}.png")
        results[name] = val_metrics
        trained_models[name] = trained
        models_for_plots[name] = (trained, val_data)
        
        print(f"\n{name} validation metrics:")
        for metric, value in val_metrics.items():
            print(f"  {metric}: {value:.4f}")
    
    # Train and validate additional baseline models
    print("\n===== Training Deep Neural Network =====")
    dnn_model = DeepNeuralNetworkModel(input_size=input_size)
    trained_dnn_model, dnn_metrics = train_filtration_AI(dnn_model, train_data)
    dnn_val_metrics = validate_filtration_AI(trained_dnn_model, val_data, output_name='DNN.png')
    results["DeepNeuralNetwork"] = dnn_val_metrics
    trained_models["DeepNeuralNetwork"] = trained_dnn_model
    models_for_plots["DeepNeuralNetwork"] = (trained_dnn_model, val_data)
    
    print("\nDeep Neural Network validation metrics:")
    for metric, value in dnn_val_metrics.items():
        print(f"  {metric}: {value:.4f}")
    
    print("\n===== Training Logistic Regression =====")
    lr_model = LogisticRegressionModel()
    trained_lr_model, lr_metrics = train_filtration_AI(lr_model, train_data)
    lr_val_metrics = validate_filtration_AI(trained_lr_model, val_data, output_name='LR.png')
    results["LogisticRegression"] = lr_val_metrics
    trained_models["LogisticRegression"] = trained_lr_model
    models_for_plots["LogisticRegression"] = (trained_lr_model, val_data)
    
    print("\nLogistic Regression validation metrics:")
    for metric, value in lr_val_metrics.items():
        print(f"  {metric}: {value:.4f}")
    
    # Generate combined plots
    print("\n===== Generating combined plots =====")
    output_dir = 'tests_outs'
    plot_combined_roc_curves(models_for_plots, output_dir=output_dir)
    print(f"Combined ROC curves saved to {os.path.join(output_dir, 'combined_roc_curves.png')}")
    
    plot_combined_metrics(results, output_dir=output_dir)
    print(f"Combined metrics plot saved to {os.path.join(output_dir, 'combined_metrics.png')}")
    
    # Select best model based on F1 score
    best_model_name = max(results.keys(), key=lambda m: results[m]["f1"])
    best_model = trained_models[best_model_name]
    
    print(f"\n{'='*60}")
    print(f"Best model: {best_model_name} (F1: {results[best_model_name]['f1']:.4f})")
    print(f"{'='*60}")
    
    # Architecture search for GAIL if it's the best model
    if best_model_name == "GAIL":
        print("\n" + "="*60)
        print("GAIL is the best model. Performing architecture search...")
        print("="*60)
        
        # Define GAIL architecture variations
        gail_variations = {
            "GAIL_Deep": lambda n: TorchClassifier(GAILDeep(n), weight_pos=5),
            "GAIL_Wide": lambda n: TorchClassifier(GAILWide(n), weight_pos=5),
            "GAIL_Narrow": lambda n: TorchClassifier(GAILNarrow(n), weight_pos=5),
            "GAIL_Dropout": lambda n: TorchClassifier(GAILWithDropout(n), weight_pos=5),
            "GAIL_Custom1": lambda n: TorchClassifier(GAILDiscriminator(n, hidden1=384, hidden2=192), weight_pos=5),
            "GAIL_Custom2": lambda n: TorchClassifier(GAILDiscriminator(n, hidden1=192, hidden2=96), weight_pos=5),
        }
        
        gail_search_results = {}
        gail_trained_models = {}
        
        for variant_name, constructor in gail_variations.items():
            print(f"\n===== Training {variant_name} =====")
            variant_model = constructor(input_size)
            
            # Train with learning curve tracking
            X_train = train_data.drop(columns=['type'])
            y_train = train_data['type']
            variant_model.train(X_train, y_train, epochs=100, batch_size=32, 
                              val_data=val_data, track_curves=True)
            
            # Validate
            variant_val_metrics = validate_filtration_AI(
                variant_model, val_data, output_name=f"{variant_name}.png"
            )
            gail_search_results[variant_name] = variant_val_metrics
            gail_trained_models[variant_name] = variant_model
            
            print(f"\n{variant_name} validation metrics:")
            for metric, value in variant_val_metrics.items():
                print(f"  {metric}: {value:.4f}")
        
        # Compare with original GAIL
        gail_search_results["GAIL_Original"] = results["GAIL"]
        gail_trained_models["GAIL_Original"] = trained_models["GAIL"]
        
        # Find best GAIL variant
        best_gail_variant = max(gail_search_results.keys(), 
                               key=lambda m: gail_search_results[m]["f1"])
        best_gail_model = gail_trained_models[best_gail_variant]
        
        print(f"\n{'='*60}")
        print(f"Best GAIL variant: {best_gail_variant} (F1: {gail_search_results[best_gail_variant]['f1']:.4f})")
        print(f"{'='*60}")
        
        # Plot learning curves for best GAIL variant
        if hasattr(best_gail_model, 'train_losses') and len(best_gail_model.train_losses) > 0:
            print("\nGenerating learning curves for best GAIL variant...")
            plot_learning_curves(best_gail_model, output_dir=output_dir, 
                               output_name=f'learning_curves_{best_gail_variant}.png')
            print(f"Learning curves saved to {os.path.join(output_dir, f'learning_curves_{best_gail_variant}.png')}")
        
        # Update best model if variant is better
        if gail_search_results[best_gail_variant]["f1"] > results["GAIL"]["f1"]:
            print(f"\nBest GAIL variant ({best_gail_variant}) outperforms original GAIL!")
            best_model = best_gail_model
            best_model_name = best_gail_variant
        else:
            print(f"\nOriginal GAIL remains the best.")
    
    # Plot learning curves for the best model (if it has learning curve data)
    if hasattr(best_model, 'train_losses') and len(best_model.train_losses) > 0:
        print("\nGenerating learning curves for best model...")
        plot_learning_curves(best_model, output_dir=output_dir, 
                           output_name='learning_curves_best_model.png')
        print(f"Learning curves saved to {os.path.join(output_dir, 'learning_curves_best_model.png')}")
    
    # Apply best model to test set
    print("\nApplying best model to test set...")
    test_predictions = apply_filtration_AI(best_model, test_data)
    
    # Save predictions
    os.makedirs(output_dir, exist_ok=True)
    test_predictions.to_csv(os.path.join(output_dir, 'test_predictions.csv'), index=False)
    print(f"\nPredictions saved to {os.path.join(output_dir, 'test_predictions.csv')}")

if __name__ == '__main__':
    main() 