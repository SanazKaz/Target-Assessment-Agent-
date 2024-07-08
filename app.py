from flask import Flask, request, jsonify
from flask_cors import CORS
from agentic import process_message

app = Flask(__name__)
CORS(app, resources={r"/chat": {"origins": "*"}})

@app.route('/chat', methods=['POST'])
def chat():
    try:
        message = request.json['message']
        conversation = process_message(message)
        return jsonify({'conversation': conversation})
    except Exception as e:
        print(f"Error processing message: {str(e)}")
        return jsonify({'error': 'An error occurred while processing your message.'}), 500

if __name__ == '__main__':
    app.run(debug=True)