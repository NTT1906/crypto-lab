import subprocess
from openai import OpenAI

# 🔑 Replace with your OpenRouter API key
client = OpenAI(
    base_url="https://openrouter.ai/api/v1",
    api_key="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",
)

def get_git_diff():
    """Get staged git diff."""
    try:
        diff = subprocess.check_output(
            ["git", "diff"], universal_newlines=True
        )
        if not diff.strip():
            print("No staged changes found.")
            exit(0)
        return diff
    except subprocess.CalledProcessError as e:
        print("Error getting git diff:", e)
        exit(1)

def generate_commit_message(diff_text):
    """Generate commit message using a free model on OpenRouter."""
    prompt = f"""
    You are a helpful assistant that writes clear and concise git commit messages.
    Based on the following diff, write a short, single-line commit message using 
    the conventional commit style (e.g., feat:, fix:, docs:, refactor:, style:, chore:).

    Git diff:
    {diff_text}
    """

    completion = client.chat.completions.create(
        model="google/gemini-2.5-flash-image",
        messages=[
            {"role": "system", "content": "You are an expert in git commit messages."},
            {"role": "user", "content": prompt},
        ],
        max_tokens=100,
        temperature=0.5,
    )

    return completion.choices[0].message.content.strip()

if __name__ == "__main__":
    diff = get_git_diff()
    commit_message = generate_commit_message(diff)
    print("\n💬 Suggested commit message:\n")
    print(commit_message)
    print("\nRun this to commit:")
    print(f'git commit -m "{commit_message}"')
